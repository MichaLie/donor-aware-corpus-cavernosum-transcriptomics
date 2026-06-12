#!/usr/bin/env python3
"""
P4 — Independent replication in GSE259348 (one DMED/diabetic-ED donor, GSM8114573),
re-fetched from GEO. Annotate pericytes ourselves, build pericyte pseudobulk, then
combine with GSE206528 pericyte pseudobulks and run a donor-level
DESeq2 ~dataset+condition; report the 8 non-stress core genes.

QC uses min_genes>=500, pct_mt<20, and scrublet. Cell types assigned
per Leiden cluster by canonical marker scores. "Pericyte/mural" = RGS5/PDGFRB-high mural.
"""
import os, numpy as np, pandas as pd, scanpy as sc, anndata as ad
from pathlib import Path
from scipy.io import mmread
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

ROOT=(os.environ.get("PROJECT_ROOT") or str(Path(__file__).resolve().parents[2]))
RAW=f"{ROOT}/data/raw/GSE259348"
OUTD=f"{ROOT}/data/processed/GSE259348"; os.makedirs(OUTD, exist_ok=True)
OUT=f"{ROOT}/results/tables/repro"; os.makedirs(OUT, exist_ok=True)
sc.settings.verbosity=1

MARKERS={
 "Pericyte":["RGS5","PDGFRB","NOTCH3","KCNJ8","HIGD1B"],
 "SMC":["ACTA2","MYH11","TAGLN","CNN1","DES"],
 "Endothelial":["PECAM1","VWF","CLDN5","CDH5"],
 "Fibroblast":["COL1A1","DCN","LUM","PDGFRA"],
 "Immune":["PTPRC","CD68","LYZ","CD3E"],
 "Epithelial":["EPCAM","KRT8","KRT18"],
}
CORE8=["IL6","ADAMTS4","PLA2G2A","CYR61","PI15","SH3BP5","ERAP2","RGS16"]

# ---- load 10x ----
print("loading 10x mtx ...")
M=mmread(f"{RAW}/GSM8114573_matrix.mtx.gz").T.tocsr()        # cells x genes
feat=pd.read_csv(f"{RAW}/GSM8114573_features.tsv.gz", sep="\t", header=None)
bc=pd.read_csv(f"{RAW}/GSM8114573_barcodes.tsv.gz", header=None)[0].values
a=ad.AnnData(M, obs=pd.DataFrame(index=bc), var=pd.DataFrame(index=feat[1].values))
a.var_names_make_unique()
a.layers["counts"]=a.X.copy()
print("raw:", a.shape)

# ---- QC ----
a.var["mt"]=a.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None)
a=a[(a.obs.n_genes_by_counts>=500) & (a.obs.pct_counts_mt<20)].copy()
print("after QC (min_genes>=500, pct_mt<20):", a.shape)
try:
    sc.pp.scrublet(a, random_state=0); a=a[~a.obs.predicted_doublet].copy()
    print("after scrublet:", a.shape)
except Exception as e:
    print("scrublet skipped:", e)

# ---- normalize / cluster ----
sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
a.raw=a
sc.pp.highly_variable_genes(a, n_top_genes=2500);
av=a[:, a.var.highly_variable].copy()
sc.pp.scale(av, max_value=10); sc.tl.pca(av, n_comps=50, random_state=0)
sc.pp.neighbors(av, n_neighbors=15, random_state=0); sc.tl.leiden(av, resolution=0.5, random_state=0, flavor="igraph", n_iterations=2, directed=False)
a.obs["leiden"]=av.obs["leiden"].values

# ---- annotate clusters by marker score ----
for ct,gs in MARKERS.items():
    sc.tl.score_genes(a, [g for g in gs if g in a.var_names], score_name=f"sc_{ct}", random_state=0)
sccols=[f"sc_{ct}" for ct in MARKERS]
clmean=a.obs.groupby("leiden", observed=True)[sccols].mean()
cl2ct={cl: sccols[int(np.argmax(row.values))].replace("sc_","") for cl,row in clmean.iterrows()}
a.obs["cell_type"]=a.obs["leiden"].map(cl2ct).astype(str)
print("\ncluster -> cell type (by max marker score):"); print(clmean.round(2).to_string())
print("\ncell type counts:", a.obs.cell_type.value_counts().to_dict())
a.write(f"{OUTD}/GSE259348_annotated_ours.h5ad")

# ---- pericyte pseudobulk (raw counts) ----
peri = a[a.obs.cell_type=="Pericyte"]
print(f"\nGSE259348 pericyte/mural cells: {peri.n_obs}")
pb259 = np.asarray(peri.layers["counts"].sum(0)).ravel()
pb259 = pd.Series(pb259, index=a.var_names, name="GSE259348_DMED")

# ---- combine with GSE206528 pericyte pseudobulk -> DESeq2 ~dataset+condition ----
cnt=pd.read_csv(f"{ROOT}/data/processed/GSE206528/pseudobulk_mural_subtypes/GSE206528_mural_subtype_pseudobulk_counts.tsv.gz", sep="\t", index_col=0)
meta=pd.read_csv(f"{ROOT}/data/processed/GSE206528/pseudobulk_mural_subtypes/GSE206528_mural_subtype_pseudobulk_metadata.tsv", sep="\t")
periM=meta[meta.mural_subtype=="Pericyte_like"].set_index("pseudobulk_id")
cols=[c for c in cnt.columns if c in periM.index]
gmap={"diabetic_ED":"diab_ED","non_diabetic_ED":"nondiab_ED","normal":"normal"}
genes=cnt.index.intersection(pb259.index)
comb=cnt.loc[genes, cols].copy(); comb["GSE259348_DMED"]=pb259.reindex(genes).fillna(0).astype(int)
md=pd.DataFrame(index=comb.columns)
md["dataset"]=["GSE206528"]*len(cols)+["GSE259348"]
md["condition"]=[gmap[periM.loc[c,"ed_subgroup"]] for c in cols]+["diab_ED"]
md["condition"]=pd.Categorical(md["condition"], categories=["normal","nondiab_ED","diab_ED"])
X=comb.T.astype(int)
keep=(X.sum(0)>=10)|X.columns.isin(CORE8); X=X.loc[:,keep]
dds=DeseqDataSet(counts=X, metadata=md, design="~dataset+condition", refit_cooks=False, quiet=True); dds.deseq2()
def ctr(a_,b_):
    st=DeseqStats(dds, contrast=["condition",a_,b_], quiet=True); st.summary()
    return st.results_df.loc[CORE8,["log2FoldChange","pvalue","padj"]]
dn=ctr("diab_ED","normal"); dne=ctr("diab_ED","nondiab_ED")
rep=pd.DataFrame({"lfc_DvsNormal":dn.log2FoldChange,"padj_DvsNormal":dn.padj,
                  "lfc_DvsNonDiabED":dne.log2FoldChange,"p_DvsNonDiabED":dne.pvalue}).loc[CORE8]
print("\n===== GSE206528 + re-fetched GSE259348 donor-level ~dataset+condition =====")
pd.set_option("display.float_format",lambda x:f"{x:.3g}")
print(rep.to_string())
rep.to_csv(f"{OUT}/04_replication_combined_DE.csv")
print(f"\nsaved -> {OUT}/04_replication_combined_DE.csv")

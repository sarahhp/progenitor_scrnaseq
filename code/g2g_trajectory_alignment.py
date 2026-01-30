import os
import re
import anndata
import scanpy as sc
import numpy as np
import seaborn as sb
import warnings
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

from genes2genes import Main
from genes2genes import ClusterUtils
from genes2genes import TimeSeriesPreprocessor
from genes2genes import PathwayAnalyser
from genes2genes import VisualUtils 


from optbinning import ContinuousOptimalBinning

output_file = snakemake.output[0]
odir = os.path.dirname(output_file)
#odir = "output/beige_vs_white/bvw_fastmnn"

path = re.split(r"[/.]",odir )
if path[0] in ["projects","home"]:
    wdir = path.index("progenitor_scrnaseq")
    module = path(wdir+1)
    dataname = path(wdir+2)
else:
    module = path[0]
    dataname = path[1]
fig_path = os.path.join("analysis", module, "figures", dataname +"_g2g"

pal = {'day0':"#B3B3B3", 'day1':"#85C2EA", 'day3':"#1F78B4"}
prog_cols = {"p0":"#DA1819","p1":"#691A93","p2":"#EBB400","p3":"#434D51",
            "p4":"#00AED1","p5":"#3AAD00","p6":"#FF7F00",
            'day0':"#B3B3B3", 'day1':"#85C2EA", 'day3':"#1F78B4"}
clusters = dict(zip([str(i) for i in range(12)],
            ["#5b859e", "#1e395f" ,"#75884b", "#1e5a46", "#df8d71", "#af4f2f" ,
            "#d48f90", "#732f30", "#ab84a5", "#59385c", "#d8b847", "#b38711"]))
bvw = ["#1F78B4","#FF7F00"]
subset_data = True
import random 
random.seed(10)

if subset_data:
    adata_ref = anndata.read_h5ad(odir + '/10%_white_monocle_pseudotime_seurat.h5ad') # Reference dataset
    adata_query = anndata.read_h5ad(odir + '/10%_beige_monocle_pseudotime_seurat.h5ad') # Query datase
elif not subset_data:
    adata_ref = anndata.read_h5ad(odir + '/white_monocle_pseudotime_seurat.h5ad') # Reference dataset
    adata_query = anndata.read_h5ad(odir + '/beige_monocle_pseudotime_seurat.h5ad') # Query datase

adata_ref.obs["exp.time"] = adata_ref.obs.time
adata_ref.obs["time"] = adata_ref.obs.monocle_pseudotime
adata_query.obs["exp.time"] = adata_query.obs.time
adata_query.obs["time"] = adata_query.obs.monocle_pseudotime

sc.pp.log1p(adata_ref)
sc.pp.highly_variable_genes(adata_ref, n_top_genes=2000)

sc.pp.log1p(adata_query)
sc.pp.highly_variable_genes(adata_query, n_top_genes=2000)
gene_list =adata_ref.var.highly_variable | adata_query.var.highly_variable
gene_list = gene_list.index[gene_list].tolist()
print(len(gene_list), "genes such as:")
print(gene_list[0:10])


x = np.asarray(adata_ref.obs.time)
optw = ContinuousOptimalBinning(name='pseudotime', dtype="numerical")
optw.fit(x, x)
print(len(optw.splits))

x = np.asarray(adata_query.obs.time)
optb = ContinuousOptimalBinning(name='pseudotime', dtype="numerical")
optb.fit(x, x)
print(len(optb.splits))

n_bins=15#int(np.mean([len(optw.splits),len(optb.splits)]))
print("Using bins", n_bins)

VisualUtils.plot_celltype_barplot(adata_ref, n_bins, "exp.time", pal)
plt.savefig(fig_path + "/exp.time_barplot.png", format="png", dpi=300, bbox_inches="tight")
plt.close()   

VisualUtils.plot_celltype_barplot(adata_ref, n_bins, "initial_clusters", clusters)
plt.savefig(fig_path +"/bvw_clusters_barplot.png", format="png", dpi=300, bbox_inches="tight")
plt.close()   

VisualUtils.plot_celltype_barplot(adata_ref, n_bins, "time.progen_clusters", prog_cols)
plt.savefig(fig_path +"/progen_clusters_barplot.png", format="png", dpi=300, bbox_inches="tight")
plt.close()   


### Run G2G alignment 
aligner = Main.RefQueryAligner(adata_ref, adata_query, gene_list, n_bins) #
aligner.align_all_pairs()

p = aligner.get_aggregate_alignment()
print(p)
plt.savefig(fig_path +"/average_heatmap.png", format="png", dpi=300, bbox_inches="tight")
plt.close()   

### Get genes
df = aligner.get_stat_df() # ordered genes according to alignment similarity statistics 
print(df.head())

diff_genes = df.iloc[0:10,0].tolist()
# for gene in diff_genes:
#     print("investigating" ,gene)
#     VisualUtils.show_gene_alignment(gene, aligner, adata_ref, adata_query, "initial_clusters", clusters)
#     plt.savefig("{}/{}_alignment_bins.png".format(fig_path, gene), format="png", dpi=300, bbox_inches="tight")
#     plt.close() 
#     
#     VisualUtils.plotTimeSeries(gene, aligner, plot_cells=True, colr=bvw)
#     plt.savefig("{}/{}_alignment.png".format(fig_path, gene), format="png", dpi=300, bbox_inches="tight")
#     plt.close() 

import pickle
if subset_data:
    with open(odir + '/g2g_aligner_10%.pkl', 'wb') as file:
        pickle.dump(aligner, file)
    df.to_csv(odir + "/10%_g2g_alignment_genes.tsv", sep='\t', index=False)
elif not subset_data:
    with open(odir + '/g2g_aligner.pkl', 'wb') as file:
        pickle.dump(aligner, file)
    df.to_csv(odir + "/g2g_alignment_genes.tsv", sep='\t', index=False)





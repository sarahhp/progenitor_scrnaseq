import os
import pickle
import anndata
import numpy as np
import seaborn as sb
import pandas as pd
import warnings
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

from genes2genes import Main
from genes2genes import ClusterUtils
from genes2genes import TimeSeriesPreprocessor
from genes2genes import PathwayAnalyser
from genes2genes import VisualUtils 

## Setup 
module = "beige_vs_white"
dataname = "bvw_fastmnn"
odir = "output/beige_vs_white/bvw_fastmnn"
fig_path = os.path.join("analysis", module, "figures", dataname +"_g2g")
subset_data = True
load_data =True

bvw = ["#1F78B4","#FF7F00"]
clusters = dict(zip([str(i) for i in range(12)],
            ["#5b859e", "#1e395f" ,"#75884b", "#1e5a46", "#df8d71", "#af4f2f" ,
            "#d48f90", "#732f30", "#ab84a5", "#59385c", "#d8b847", "#b38711"]))

## Open previous alignment results
os.makedirs(fig_path, exist_ok=True)
with open(odir + '/g2g_aligner.pkl', 'rb') as file:
    aligner = pickle.load(file)
    
## Plot a gene
genes=[
    #R1 "PEMT","PLIN1","ICAM1","PM20D1","MGP"
    #"FABP4","FABP3","FASN","ACLY","LPL","PCK1",
    #"BRD4","SLC3A2",
    #"VIM","ACTA2","TAGLN","PHLDA1",
    "CD36","CXCL8"
    ]
for gene in genes:
    if gene not in aligner.results_map.keys():
        print(gene,"not found in alignment result. Is gene not in HVG list or mispelled? Skipping")
    else:
        print("Plotting " + gene)
        VisualUtils.plotTimeSeries(gene, aligner, plot_cells=True)
        plt.savefig("{}/{}_alignment.png".format(fig_path, gene), format="png", dpi=300, bbox_inches="tight")
        plt.close() 
        
print("Plotting without cells" , gene)
VisualUtils.plotTimeSeries(gene, aligner, plot_cells=False)
plt.savefig("{}/{}_alignment.png".format(fig_path, gene), format="png", dpi=300, bbox_inches="tight")
plt.close() 

## Open and process anndata files
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
#sc.pp.log1p(adata_ref)
#sc.pp.log1p(adata_query)

## Barplot
print("Silently plotting basic barplot; which adds bin ids to adata")
VisualUtils.plot_celltype_barplot(adata_ref, 15, "initial_clusters", clusters)
VisualUtils.plot_celltype_barplot(adata_query, 15, "initial_clusters", clusters)
print("Trying barplot" , gene)
VisualUtils.visualize_gene_alignment(aligner.results_map[gene], adata_ref, adata_query, 
"initial_clusters", cmap=clusters)
#VisualUtils.show_gene_alignment(gene, aligner, adata_ref, adata_query, "initial_clusters", clusters)
plt.savefig("{}/{}_alignment_bins.png".format(fig_path, gene), format="png", dpi=300, bbox_inches="tight")
plt.close()


adata_ref.obs['bin_ids'].to_csv(os.path.join(odir, "g2g_white_bin_ids.tsv"), sep='\t')
adata_query.obs['bin_ids'].to_csv(os.path.join(odir, "g2g_beige_bin_ids.tsv"), sep='\t')



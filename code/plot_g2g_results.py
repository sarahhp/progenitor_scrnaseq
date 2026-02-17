import os
import pickle
import anndata
import numpy as np
import seaborn as sb
import pandas as pd
import warnings
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
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
    "ACTA2","TAGLN","PHLDA1","MGP","PLIN1","PCK1",
    "PEMT","ICAM1","PM20D1",
    "FABP4","FABP3","FASN","ACLY","LPL",
    "BRD4","SLC3A2",
    "VIM",
    "CD36","CXCL8",
    #diff genes
    "ENPP5", "ITPKA", "FPR2",
    "IL23A",   "AGXT" , "GALNT14" , "TSPAN19" , "CA4" , "MMP7",  "SSTR2",   "HIST1H2BI",
    "DIO2"
    ]
for gene in genes:
    if gene not in aligner.results_map.keys():
        print(gene,"not found in alignment result. Is gene not in HVG list or mispelled? Skipping")
    else:
        print("Plotting " + gene)
        VisualUtils.plotTimeSeries(gene, aligner, plot_cells=True)
        #Keep only left-most figure (trajectory)
        axes = plt.gcf().get_axes()
        axes[1].remove()
        axes[2].remove()
        
        #Re-add x axis
        left_ax = axes[0]
        left_ax.set_axis_on()
        left_ax.set_xlabel("Pseudotime (Monocle)")
        left_ax.set_ylabel("TSS expression")
        
        #change colours
        colors = ["#1F78B4","#FF7F00"]
        
        for line, col in zip(left_ax.lines, colors):
            line.set_color(col)
        
        scatters = [c for c in left_ax.collections if hasattr(c, "get_offsets")]
        for coll, col in zip(scatters, colors):
            alpha = coll.get_alpha() or 1.0 #keep previous alpha
            coll.set_facecolor(mcolors.to_rgba(col, alpha))
        for device in ["png","pdf"]:
            plt.savefig("{}/{}_alignment.{}".format(fig_path, gene, device), format=device, dpi=300, bbox_inches="tight")
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

print("Silently plotting basic barplot; which adds bin ids to adata")
VisualUtils.plot_celltype_barplot(adata_ref, 15, "initial_clusters", clusters)
VisualUtils.plot_celltype_barplot(adata_query, 15, "initial_clusters", clusters)

## Barplot
for gene in genes:
    if gene not in aligner.results_map.keys():
        print(gene,"not found in alignment result. Is gene not in HVG list or mispelled? Skipping")
    else:
        print("Plotting barplot" , gene)
        VisualUtils.visualize_gene_alignment(aligner.results_map[gene], adata_ref, adata_query, 
        "initial_clusters", cmap=clusters)
        plt.savefig("{}/{}_alignment_bins.png".format(fig_path, gene), format="png", dpi=300, bbox_inches="tight")
        plt.close()


adata_ref.obs['bin_ids'].to_csv(os.path.join(odir, "g2g_white_bin_ids.tsv"), sep='\t')
adata_query.obs['bin_ids'].to_csv(os.path.join(odir, "g2g_beige_bin_ids.tsv"), sep='\t')



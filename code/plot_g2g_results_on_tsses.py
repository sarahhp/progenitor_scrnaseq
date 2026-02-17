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
dataname = "tss_fastmnn"
odir = "output/beige_vs_white/tss_fastmnn"
fig_path = os.path.join("analysis", module, "figures", dataname +"_g2g")
subset_data = True
load_data =True

bvw = ["#1F78B4","#FF7F00"]
clusters = dict(zip([str(i) for i in range(12)],
            ["#417BA1", "#ff1493", "#ffff2e", "#345634" ,"#8b0000",
            "#ff6700", "#93C0D5" ,"#8b4513", "#9248A7", "#1c8859",
            "#474747", "#8fbc8f", "#d2b48c"]))

## Open previous alignment results
os.makedirs(fig_path, exist_ok=True)
with open(odir + '/g2g_aligner.pkl', 'rb') as file:
    aligner = pickle.load(file)
    
## Plot a TSS
genes=[
# "P:r2@BRD4", "P:r1@PROX1", "P:r1@ENPP5",
 #"P:r3@EPB41L3", "P:r1@PPP1R1A", "P:r2@GREB1L",
 #"P:r1@CLDN7", "P:r1@PCSK2", "P:r1@KCNB1",
 #"P:r1@TSPAN19", "P:r3@PEMT", "P:r1@PEMT",
 "P:r1@PLIN1", #"P:r1@ICAM1", "P:r1@PM20D1",
 "P:r1@MGP", #"P:r1@FABP4", "P:r1@FABP3",
 "P:r1@LPL", #"P:r1@PCK1", "P:r1@ACTA2",
 "P:r1@TAGLN", "P:r1@PHLDA1",# "P:r3@CD36",
 #"P:r6@CD36", "P:r1@CD36", "P:r2@CD36",
 #"P:r1@CXCL8", "P:r1@ITPKA", "P:r1@FPR2",
 #"P:r1@AGXT", "P:r1@CA4", "P:r1@MMP7",
 #"P:r1@HIST1H2BI", "P:r1@DIO2"
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
        colors = ["#FF7F00", "#1F78B4"]
        
        for line, col in zip(left_ax.lines, colors):
            line.set_color(col)
        
        scatters = [c for c in left_ax.collections if hasattr(c, "get_offsets")]
        for coll, col in zip(scatters, colors):
            alpha = coll.get_alpha() or 1.0 #keep previous alpha
            coll.set_facecolor(mcolors.to_rgba(col, alpha))

        with open("{}/{}_alignment_fig.pkl".format(fig_path, gene), 'wb') as file:
            pickle.dump(plt.gcf(), file)
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



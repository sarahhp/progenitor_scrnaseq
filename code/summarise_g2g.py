"""Compare trajectories using g2g alignment tool  """

__author__ = "Sarah Hazell Pickering (s.h.pickering@medisin.uio.no)"
__date__ = "2026-01-12"


import re
include: "run_g2g.py"

def choose_input(wildcards):
    ''' NB this syntax doesn't allow  for cross-factor comparisons e.g.
    monocle white vs scvelo beige 
    '''
    print("Choosing input files for rule g2g_comparsions")
    comp = wildcards.comp
    if "beige_vs_white" in comp:
        method = re.sub("_beige_vs_white", "", comp) #convert to python 
        if method == "monocle":
            input1 = wildcards.dir + "/white_monocle_pseudotime_seurat.h5ad"
            input2 = wildcards.dir + "/beige_monocle_pseudotime_seurat.h5ad"
            return([input1,input2])
        elif method == "scvelo":
            input1 = wildcards.dir + "/dynamical_velocity_scexperiment_white.h5ad"
            input2 = wildcards.dir + "/dynamical_velocity_scexperiment_beige.h5ad"
            return([input1,input2])
    if "monocle_vs_scvelo" in comp:
        print("Monocle vs scvelo comparison, checking dataset")
        dataset = re.sub("_monocle_vs_scvelo", "", comp) 
        if dataset == "beige":
            input1 =  wildcards.dir + "/beige_monocle_pseudotime_seurat.h5ad"
            input2 =  wildcards.dir + "/dynamical_velocity_scexperiment_beige.h5ad"
            return([input1,input2])
        elif dataset == "white":
            print("dataset = white")
            input1 = wildcards.dir + "/white_monocle_pseudotime_seurat.h5ad"
            input2 = wildcards.dir + "/dynamical_velocity_scexperiment_white.h5ad"
            return([input1,input2])
    else:
        print("File pattern not recognised. Wildcards", wildcards)

rule all:
    input:
        expand("{dir}/g2g_alignments/{comp}/{comp}_g2g.log",
            dir = "output/beige_vs_white/bvw_fastmnn",
            comp = "monocle_beige_vs_white")

#e.g. monocle input
# "{dir}/white_monocle_pseudotime_seurat.h5ad",
#         "{dir}/beige_monocle_pseudotime_seurat.h5ad"


rule g2g_comparisons:
    input:
        anndata = choose_input
    conda: "g2g_py3.9"
    log: "{dir}/g2g_alignments/{comp}/{comp}_g2g.log"
    wildcard_constraints:
        comp="[^/]+" #cannot contain "/"
    output:
        #"{dir}/g2g_alignments/{comp}/g2g_aligner.pkl",
        "{dir}/g2g_alignments/{comp}/g2g_alignment_genes.tsv"
    script:
        "g2g_align_any_trajectory.py"
        
        
def choose_comparison(wildcards):
    return wildcards.dir

rule g2g_across:
    input:
        anndata = choose_comparison,
        script = "code/g2g_align_any_trajectory.py"
    conda: "g2g_py3.9"
    log: "{dir}/bvw_g2g_alignments/fastmnn_vs_harmony/{comp1}/{comp2}_g2g.log"
    wildcard_constraints:
        comp="[^/]+" #cannot contain "/"
    output:
       # "{dir}/bvw_g2g_alignments/fastmnn_vs_harmony/{comp1}/{comp2}_g2g_aligner.pkl",
        "{dir}/bvw_g2g_alignments/fastmnn_vs_harmony/{comp1}/{comp2}_g2g_alignment_genes.tsv"
    shell:
        "python {input.script} > {log}"
        
rule g2g_heatmap:
    input:
        within_method = expand("{odir}/bvw_{method}/g2g_alignments/{comp}/{comp}_g2g.log",
                odir = "output/beige_vs_white",
                method= ["fastmnn","harmony"],
                comp=["monocle_beige_vs_white","scvelo_beige_vs_white", 
                        "beige_monocle_vs_scvelo","white_monocle_vs_scvelo"]),
        across_methods = expand("{odir}/bvw_g2g_alignments/fastmnn_vs_harmony/{comp1}/{comp2}_g2g.log",
                odir = "output/beige_vs_white",
                comp1=["monocle_beige","scvelo_beige","monocle_white", "scvelo_white"],
                comp2=["monocle_beige","scvelo_beige","monocle_white", "scvelo_white"])
    output:
        expand("{odir}/bvw_g2g_alignments/heatmap_summary.pdf",
                odir = "output/beige_vs_white")
    
                
                
                
                
        

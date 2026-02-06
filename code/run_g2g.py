"""   """

__author__ = "Sarah Hazell Pickering (s.h.pickering@medisin.uio.no)"
__date__ = "2025-08-20"

import os

rule test:
    input:
        expand("{dir}/{file}.h5ad",
                dir = "output/beige_vs_white/bvw_fastmnn",
                file = "10%_white_monocle_pseudotime_seurat")

rule seurat_to_anndata:
    ''' Version 2 loads the RStudio module and checks if it needs
    to create the reticulate conda env sceasy.
    '''
    input:
        rdata = "{dir}/{module}/{dataname}/{file}.rds",
        rscript = "code/convert_seurat_to_anndata.R"
    params:
        condaenv = "sceasy",
        MODULE = "RStudio/2023.12.1+402-1-R-4.2.1"
    output:
        anndata = "{dir}/{module}/{dataname}/{file}.h5ad",
        subset = "{dir}/{module}/{dataname}/10%_{file}.h5ad",
    shell:
        "if ! conda env list | grep -q \"^{params.condaenv}\"; then "
            "echo 'Creating conda env: sceasy' ; "
            "mamba env create --file envs/sceasy.yml -y; fi ; "
        "if ! module -t list 2>&1 | grep -q '^${params.MODULE}$'; then "
            "echo 'Loading module: {params.MODULE}' ; "
            "module load {params.MODULE}; fi ; "
        "Rscript  {input.rscript} {input.rdata} {output.anndata}"  
        

rule g2g_alignment:
    input:
        "{dir}/{is_sbs}white_monocle_pseudotime_seurat.h5ad",
        "{dir}/{is_sbs}beige_monocle_pseudotime_seurat.h5ad"
    conda: "g2g_py3.9"
    wildcard_constraints:
        is_sbs="(10%_)*"
    output:
        "{dir}/{is_sbs}g2g_aligner.pkl",
        "{dir}/{is_sbs}g2g_alignment_genes.tsv"
    script:
        "g2g_trajectory_alignment.py"


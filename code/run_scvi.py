"""   """

__author__ = "Sarah Hazell Pickering (s.h.pickering@medisin.uio.no)"
__date__ = "2025-08-20"

import os

SCVI_DIR = "output/progenitors/scvi_integration"
#INDIR = #output/progenitors/initial"

def choose_indataname(wildcards):
    if wildcards.module == "progenitors":
        indataname = "initial"
    elif wildcards.module == "beige_vs_white":
        if wildcards.dataname.startswith("tss"):
            indataname = "tss_bvw_initial"
        elif wildcards.dataname.startswith("bvw"):
            indataname = "downsample"
    elif wildcards.module == "adipogenesis":
        indataname = "white_only_downsample"

    return(os.path.join(wildcards.dir, wildcards.module,indataname,"complete_analysis.rds"))
    

rule all:
    input:
        expand("{dir}/{file}.h5ad",
                dir = SCVI_DIR,
                file = "10%_10%_complete_analysis.scviintegrated")

rule seurat_to_anndata:
    ''' Version 2 loads the RStudio module and checks if it needs
    to create the reticulate conda env sceasy.
    '''
    input:
        rdata = choose_indataname,
        rscript = "code/convert_seurat_to_anndata.R"
    params:
        condaenv = "sceasy",
        MODULE = "RStudio/2023.12.1+402-1-R-4.2.1"
    output:
        anndata = "{dir}/{module}/{dataname}/complete_analysis.h5ad",
        subset = "{dir}/{module}/{dataname}/10%_complete_analysis.h5ad",
    shell:
        "if ! conda env list | grep -q \"^{params.condaenv}\"; then "
            "echo 'Creating conda env: sceasy' ; "
            "mamba env create --file envs/sceasy.yml -y; fi ; "
        "if ! module -t list 2>&1 | grep -q '^${params.MODULE}$'; then "
            "echo 'Loading module: {params.MODULE}' ; "
            "module load {params.MODULE}; fi ; "
        "Rscript  {input.rscript} {input.rdata} {output.anndata}"  
        

rule scvi:
    input:
        "{dir}/{file}.h5ad"
    conda: "../envs/scvi.yml"
    output:
        "{dir}/{file}.scviintegrated.h5ad"
    script:
        "scvi_integration.py"


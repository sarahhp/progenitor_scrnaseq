"""   """

__author__ = "Sarah Hazell Pickering (s.h.pickering@medisin.uio.no)"
__date__ = "2025-08-20"

OUT = "output/progenitors"
INDATA = "initial"
DATANAME = "scvi_integration"

rule all:
    input:
        expand("{dir}/{dataname}/{file}.h5ad",
                dir = OUT,
                dataname=DATANAME,
                file = "10%_10%_complete_analysis.scviintegrated")

rule seurat_to_anndata:
    ''' Version 2 loads the RStudio module and checks if it needs
    to create the reticulate conda env sceasy.
    '''
    input:
        rdata = expand("{dir}/{indata}/{{file}}.rds",
                dir = OUT,
                indata = INDATA),
        rscript = "code/convert_seurat_to_anndata.R"
    params:
        condaenv = "sceasy",
        MODULE = "RStudio/2023.12.1+402-1-R-4.2.1"
    output:
        anndata = expand("{dir}/{dataname}/{{file}}.h5ad",
                dir = OUT,
                dataname=DATANAME),
        subset = expand("{dir}/{dataname}/10%_{{file}}.h5ad",
                dir = OUT,
                dataname=DATANAME)
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
        expand("{dir}/{dataname}/{{file}}.h5ad",
                dir = OUT,
                dataname=DATANAME)
    conda: "../envs/scvi.yml"
    output:
        expand("{dir}/{dataname}/{{file}}.scviintegrated.h5ad",
                dir = OUT,
                dataname=DATANAME)
    script:
        "scvi_integration.py"


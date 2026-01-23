"""   """

__author__ = "Sarah Hazell Pickering (s.h.pickering@medisin.uio.no)"
__date__ = "2025-08-20"

SCVI_DIR = config["scvi_dir"]#"output/progenitors/scvi_integration"
INDIR = config["indir"]#output/progenitors/initial"

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
        rdata = expand("{indir}/{{file}}.rds",
                indir = INDIR),
        rscript = "code/convert_seurat_to_anndata.R"
    params:
        condaenv = "sceasy",
        MODULE = "RStudio/2023.12.1+402-1-R-4.2.1"
    output:
        anndata = expand("{dir}/{{file}}.h5ad",
                dir = SCVI_DIR),
        subset = expand("{dir}/10%_{{file}}.h5ad",
                dir = SCVI_DIR)
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
        expand("{dir}/{{file}}.h5ad",
                dir = SCVI_DIR)
    conda: "../envs/scvi.yml"
    output:
        expand("{dir}/{{file}}.scviintegrated.h5ad",
                dir = SCVI_DIR)
    script:
        "scvi_integration.py"


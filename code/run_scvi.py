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
    input:
        expand("{dir}/{indata}/{{file}}.rds",
                dir = OUT,
                indata = INDATA)

    conda: "../envs/sceasy.yml"
    output:
        expand("{dir}/{dataname}/{{file}}.h5ad",
                dir = OUT,
                dataname=DATANAME),
        expand("{dir}/{dataname}/10%_{{file}}.h5ad",
                dir = OUT,
                dataname=DATANAME)
    script:
        "convert_rds_to_AnnData.R"

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


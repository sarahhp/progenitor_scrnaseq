library(sceasy)
library(BPCells)
library(Seurat)
library(reticulate)
library(here); i_am("code/rscript_sceasy.R")
use_condaenv('sceasy')

input_file = "output/progenitors/initial/complete_analysis.rds" #snakemake@input[[1]] 
output_file = "output/progenitors/scvi_integration/complete_analysis.h5ad" #snakemake@output[[1]]

print("Loading seurat object")
seurat_object = readRDS(here(input_file))
DefaultAssay(seurat_object) = "RNA"

print("Converting to in-memory")
seurat_object[["RNA"]]$scale.data = as(seurat_object[["RNA"]]$scale.data, Class="matrix")

print("Converting to v4/3 assay")
seurat_object[["RNA3"]] = as(seurat_object[["RNA"]],  Class="Assay")

print("Reassigning")
DefaultAssay(seurat_object) = "RNA3"
seurat_object[["RNA"]] = NULL
seurat_object = RenameAssays(seurat_object, RNA3="RNA")
cat("Converting full object and saving as")
cat(output_file)
cat("\n")
sceasy::convertFormat(seurat_object, from="seurat", to="anndata",
                       outFile=here(output_file), main_layer="counts")

sbs_file = here(dirname(output_file), paste0("10%_", basename(output_file)))

Idents(seurat_object) = "orig.ident"
sbs = subset(seurat_object, downsample=1000)
print("Converting... subset")
sceasy::convertFormat(sbs, from="seurat", to="anndata",
                      outFile=sbs_file,  main_layer="counts")

print("Completed!")




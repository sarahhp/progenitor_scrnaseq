import os
import tempfile

import scanpy as sc
import scvi
import seaborn as sns
import torch
#from rich import print
#from scib_metrics.benchmark import Benchmarker

input_file = snakemake.input[0]
output_file = snakemake.output[0]
odir = os.path.dirname(output_file)

print("Loading ", input_file)
adata = sc.read(input_file)
print(adata)
print(adata.layers.keys())

print("setting up scvi")
scvi.model.SCVI.setup_anndata(adata, batch_key="orig.ident")

model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")

print("Training")
model.train(batch_size=256) #slightly larger batch size to speed up learning

print("Saving the model")
model.save(os.path.join(odir, "scvi_model"), overwrite=True)

print("Get latent variables")
SCVI_LATENT_KEY = "X_scVI"
adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()

SCVI_NORMALIZED_KEY = "scvi_normalized"
adata.layers[SCVI_NORMALIZED_KEY] = model.get_normalized_expression(library_size=10e4)

print("Write")
adata.write(output_file)

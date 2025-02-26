# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.6
#   kernelspec:
#     display_name: Python 3
#     name: python3
# ---

# %% [markdown]
# # This notebook was run on a Google Colab environment with GPU acceleration enabled

# %% colab={"base_uri": "https://localhost:8080/"} id="Djfn1UCfzvbi" outputId="d9611831-f9e4-4bcc-b9a5-a4544ade2e23"
# !pip install scvi-tools scanpy anndata

# %% id="u18ngkd90ADn"
import scanpy as sc
import scvi
import anndata as ad
import numpy as np

# %% colab={"base_uri": "https://localhost:8080/"} id="_l81KWbDHXEK" outputId="e3086202-52e6-4095-fc85-03a3e1c2502f"
from google.colab import drive
drive.mount('/content/drive')

# %% [markdown]
# # Load the data

# %% id="16yvo7VwHdkt"
rna_path = "./drive/MyDrive/scRNA_all.h5ad"
adata_rna = sc.read_h5ad(rna_path)

atac_path = "./drive/MyDrive/scATAC_all_gene_activity.h5ad"
adata_atac = sc.read_h5ad(atac_path)
adata_atac.layers['lognormalised'] = adata_atac.X.copy()

adata_atac.X = adata_atac.layers['count'].copy()
adata_rna.X = adata_rna.layers['counts'].copy()

# %% [markdown]
# # Subset to only highly variable genes

# %% id="pkIapbX_JbCF"
hvg = adata_rna.var.query('highly_variable').index.tolist()
atac_genes = adata_atac.var.index.tolist()
common_genes = list(set(hvg) & set(atac_genes))

adata_rna = adata_rna[:, common_genes].copy()
adata_atac = adata_atac[:, common_genes].copy()

# %% [markdown]
# # Wrangle cell type label and combine scRNA and scATAC data

# %% id="0qU3XyOZpZYE"
# Assign batch labels
adata_rna.obs["batch"] = "scRNA"
adata_atac.obs["batch"] = "scATAC"

# Ensure cell type labels exist in scRNA-seq
adata_rna.obs["cell_type"] = adata_rna.obs["cell_type"].astype("category")
adata_rna.obs["cell_type_fine"] = adata_rna.obs["cell_type_fine"].astype("category")

# Assign "Unknown" labels to scATAC-seq
adata_atac.obs["cell_type"] = "Unknown"
adata_atac.obs["cell_type_fine"] = "Unknown"

# Concatenate both datasets
adata_combined = ad.concat([adata_rna, adata_atac], join="inner")
adata_combined.obs["batch"] = adata_combined.obs["batch"].astype("category")
adata_combined.obs["cell_type"] = adata_combined.obs["cell_type"].astype("category")
adata_combined.obs["cell_type_fine"] = adata_combined.obs["cell_type_fine"].astype("category")


# %% [markdown]
# # Set up SCVI for combined dataset

# %% colab={"base_uri": "https://localhost:8080/", "height": 223, "referenced_widgets": ["48ea7cc828724feca3bf14e3180b49fd", "e515e0905e5c40af984a2f9973dc145f", "be2138bf6b064d83b5a8372969531caf", "92e985afd9914530a0c345c8c290b2da", "9c162671ef66455f9b119064ba2bd36c", "c75d027dbb0c437697d60c78c5418ea6", "6e3754b156dc4a90a9f93b0f6d33ddfd", "58115aa9321a42069505fba0fca5e67e", "cdf166d0609c48ac9bd0418ef26a8ea1", "a839342841614d9cb4125481abc13f1e", "1b51ef3538da42c9bd7eff61afcefefa"]} id="oFyayrUGJM_f" outputId="e6854e7b-3206-4c92-d8e6-f08f392c764b"
scvi.model.SCVI.setup_anndata(adata_combined, batch_key="batch", labels_key="cell_type_fine")

# Train SCVI on both scRNA + scATAC
vae = scvi.model.SCVI(
    adata_combined,
    n_latent=30,
    n_layers=2,
    n_hidden=128,
    dropout_rate=0.2
)
vae.train(max_epochs=100, batch_size=256, plan_kwargs={"lr": 1e-3})

# %% colab={"base_uri": "https://localhost:8080/", "height": 449} id="2EGiyq1Zbhay" outputId="474b527b-65d7-41cf-d885-5c17803a9bcb"
import matplotlib.pyplot as plt

# Plot training loss
plt.plot(vae.history["elbo_train"], label="Train Loss")
plt.xlabel("Epoch")
plt.ylabel("ELBO Loss")
plt.legend()
plt.show()

# %% [markdown]
# # Initialize scANVI from trained SCVI model

# %% colab={"base_uri": "https://localhost:8080/", "height": 240, "referenced_widgets": ["a4043426c18341fc842c6e9507c85c90", "9798fc6000df40acac241e15a306ad19", "73d379792247495d973e5fccbd0fbccc", "ef42d8f39c5244838139d402e4d4fb1b", "b9597846a8534aba9d9a5e6562ccd395", "853fce1604df4efb912225b19d4141a1", "d26ee711fee4401ab575a01dacc807e7", "5d8513faf1e34ef0a2ebce40b806cf2c", "2fc89664b2ca40e09ffe51aa02c75f1d", "6c698f7f9c864b45b29b0136b271bd2e", "3297b89e0e6c4cfcb8bf7a6c05e0fef3"]} id="G7pXysWKIltp" outputId="d2416e53-4c19-43a5-9d27-8e2417aadd73"
scanvi = scvi.model.SCANVI.from_scvi_model(vae, unlabeled_category="Unknown")

# Train scANVI (semi-supervised learning)
scanvi.train(max_epochs=50, batch_size=256, plan_kwargs={"lr": 1e-3})

# %% id="YxGKIYtFrstW"
# Predict labels for scATAC-seq
adata_combined.obs["predicted_labels"] = scanvi.predict()

# Extract scATAC-seq predictions
adata_atac.obs["predicted_labels"] = adata_combined.obs["predicted_labels"].loc[adata_atac.obs_names]

# %% colab={"base_uri": "https://localhost:8080/", "height": 448} id="dCc0O1q7fR1s" outputId="db238f7f-bc0b-4405-f0e7-e35707fbe05a"
sc.pl.umap(adata_atac, color='predicted_labels')

# %% [markdown]
# # Save annnotaed scATAC-seq object and model weights

# %% id="hG8_1Rs_fXXH"
atac_path = "./drive/MyDrive/scATAC_all_gene_activity.h5ad"
adata_atac.write_h5ad(atac_path, compression='gzip')

# %% id="gmV2OexXz1X4"
import os

# Define the Google Drive save directory
save_dir = "./drive/MyDrive/scVI_model"

# Create the directory if it doesn't exist
os.makedirs(save_dir, exist_ok=True)


# %% colab={"base_uri": "https://localhost:8080/"} id="T4qQI6PZ0Q2E" outputId="1f3b5da5-0e8c-4a86-dc87-b05c6436a6dd"
vae.save(os.path.join(save_dir, "scVI_model"), overwrite=True)
print(f"Saved scVI model to {save_dir}")

# %% colab={"base_uri": "https://localhost:8080/"} id="-TzrjGgq0S7v" outputId="9f748e34-e277-430e-c975-7d0388d8d3d0"
scanvi.save(os.path.join(save_dir, "scANVI_model"), overwrite=True)
print(f"Saved scANVI model to {save_dir}")


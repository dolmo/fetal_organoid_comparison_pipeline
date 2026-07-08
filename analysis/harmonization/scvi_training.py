"""Highly-variable gene selection and scVI model training."""

import scanpy as sc
import scvi

import config


def select_hvgs(adata_fetal, adata_organoid, n_top_genes=None, flavor=None):
    """Select the union of highly-variable genes across both datasets so we
    capture drivers of variation from the fetal and organoid systems."""
    if n_top_genes is None:
        n_top_genes = config.N_TOP_HVGS
    if flavor is None:
        flavor = config.HVG_FLAVOR

    sc.pp.highly_variable_genes(adata_fetal, n_top_genes=n_top_genes, flavor=flavor)
    sc.pp.highly_variable_genes(adata_organoid, n_top_genes=n_top_genes, flavor=flavor)

    hvg_fetal = adata_fetal.var[adata_fetal.var["highly_variable"]].index
    hvg_organoid = adata_organoid.var[adata_organoid.var["highly_variable"]].index
    hvg_union = hvg_fetal.union(hvg_organoid)

    adata_fetal = adata_fetal[:, hvg_union].copy()
    adata_organoid = adata_organoid[:, hvg_union].copy()
    return adata_fetal, adata_organoid


def train_scvi(adata, batch_key=None):
    """Train an scVI model and store the batch-corrected latent representation
    in ``adata.obsm['X_scVI']``. The batch key fixes internal (per-individual)
    variation while leaving the cross-species signal intact."""
    if batch_key is None:
        batch_key = config.SCVI_BATCH_KEY

    scvi.model.SCVI.setup_anndata(adata, batch_key=batch_key)
    model = scvi.model.SCVI(
        adata,
        n_latent=config.N_LATENT,
        n_layers=config.N_LAYERS,
        dispersion=config.DISPERSION,
    )
    model.train(
        max_epochs=config.MAX_EPOCHS,
        batch_size=config.BATCH_SIZE,
        plan_kwargs={"lr": config.LEARNING_RATE},
        early_stopping=config.EARLY_STOPPING,
        early_stopping_patience=config.EARLY_STOPPING_PATIENCE,
    )
    adata.obsm["X_scVI"] = model.get_latent_representation()
    return adata

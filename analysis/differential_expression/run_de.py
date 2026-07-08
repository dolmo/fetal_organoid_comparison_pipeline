#!/usr/bin/env python3
"""
Differential expression validation stage (Phase 3).

Independently validates the RIMA / COPE findings by testing whether the ECM
("SLUG-noid", ``ecmy``) organoids up-regulate the key neurodevelopmental and
maturation genes relative to controls (``ecmn``). Runs a Wilcoxon differential
expression test on the full organoid dataset and produces a targeted dot plot.

Previously lived at workflow/diffanalysis/slugnoid_validation.py; consolidated
here so all analysis code lives under analysis/.

Run:
    python3 run_de.py
"""

import os

import scanpy as sc
import matplotlib.pyplot as plt

# ==========================================
# 1. CONFIGURATION
# ==========================================
# Pointing directly to the full organoid dataset.
ORGANOID_PATH = "/data/input/rnh027_complete_08072025.h5ad"
OUTPUT_DIR = "/data/outputs"

# Top hits from the COPE analysis plus the genes named in the abstract's
# hypothesis, grouped by biological role.
TARGET_GENES = [
    "SCN2A", "PCLO", "LRRC7",      # Synaptic & action potential machinery
    "COL4A5", "COL4A6", "LAMB2",   # ECM & structural laminin
    "SLC1A3", "KCNJ10",            # Astrocyte maturation
    "GLI3", "NOTCH1",              # Morphogens & patterning
    "SATB2", "MEF2C",              # Neuronal maturation (from abstract)
    "CARD16", "GSDMD",             # Top general COPE hits
]


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # ==========================================
    # 2. LOAD DATA
    # ==========================================
    print("Loading organoid dataset...")
    adata = sc.read_h5ad(ORGANOID_PATH)

    # Ensure the genes exist in the dataset to avoid KeyError crashes.
    genes_to_plot = [g for g in TARGET_GENES if g in adata.var_names]

    # ==========================================
    # 3. DIFFERENTIAL EXPRESSION (SLUG-noids vs Controls)
    # ==========================================
    print("Calculating differential expression (ecmy vs ecmn)...")
    sc.tl.rank_genes_groups(
        adata,
        groupby="ecm",
        groups=["ecmy"],
        reference="ecmn",
        method="wilcoxon",
    )

    de_results = sc.get.rank_genes_groups_df(adata, group="ecmy")
    csv_path = os.path.join(OUTPUT_DIR, "SLUGnoid_vs_Control_DE_Stats.csv")
    de_results.to_csv(csv_path, index=False)
    print(f"Saved Differential Expression statistics to {csv_path}")

    # ==========================================
    # 4. VISUALIZATION: DOT PLOT
    # ==========================================
    print("Generating Dot Plot for targeted maturation genes...")
    # standard_scale="var" normalizes colors from 0 to 1 per gene, making it
    # easy to see which condition has higher expression.
    sc.pl.dotplot(
        adata,
        var_names=genes_to_plot,
        groupby="ecm",
        standard_scale="var",
        title="Molecular Specification: SLUG-noids (ecmy) vs Controls (ecmn)",
        cmap="Blues",
        show=False,
    )

    plot_path = os.path.join(OUTPUT_DIR, "SLUGnoid_Targeted_Maturation_DotPlot.pdf")
    plt.savefig(plot_path, bbox_inches="tight")
    plt.close()
    print(f"Dot plot saved successfully to {plot_path}")
    print("Validation complete!")


if __name__ == "__main__":
    main()

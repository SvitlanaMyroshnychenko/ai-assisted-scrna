#!/usr/bin/env python
"""Config-driven Scanpy pipeline for remote scRNA-seq analysis.

This script is intended to be executed on the remote analysis server where the
real sample sheet and input count matrices are available. The local repository
contains templates only.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Any

import anndata as ad
import pandas as pd
import scanpy as sc
import yaml


DEFAULT_CONFIG = Path("configs/pipeline.template.yaml")
DEFAULT_TEMPLATE_SAMPLE_SHEET = Path("configs/samples.template.csv")
REQUIRED_SAMPLE_COLUMNS = ("sample_id", "path", "condition", "batch")
BIOLOGICAL_HARMONY_KEYS = {"condition", "disease", "treatment", "group", "phenotype"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a reproducible Scanpy scRNA-seq analysis pipeline."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help="Path to pipeline YAML config.",
    )
    parser.add_argument(
        "--sample-sheet",
        type=Path,
        default=None,
        help=(
            "Optional path to sample metadata CSV. If omitted, the template "
            "config uses configs/samples.template.csv; other configs use "
            "input.sample_sheet from the YAML."
        ),
    )
    return parser.parse_args()


def setup_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def read_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict):
        raise ValueError(f"Config must contain a YAML mapping: {path}")
    return config


def resolve_sample_sheet(config_path: Path, config: dict[str, Any], arg_path: Path | None) -> Path:
    if arg_path is not None:
        return arg_path

    if config_path.as_posix() == DEFAULT_CONFIG.as_posix():
        return DEFAULT_TEMPLATE_SAMPLE_SHEET

    sample_sheet = config.get("input", {}).get("sample_sheet")
    if not sample_sheet:
        raise ValueError("Missing input.sample_sheet in config.")
    return Path(sample_sheet)


def read_samples(path: Path) -> pd.DataFrame:
    samples = pd.read_csv(path)
    missing = [column for column in REQUIRED_SAMPLE_COLUMNS if column not in samples.columns]
    if missing:
        raise ValueError(f"Sample sheet is missing required columns: {missing}")
    if samples["sample_id"].duplicated().any():
        duplicated = samples.loc[samples["sample_id"].duplicated(), "sample_id"].tolist()
        raise ValueError(f"Sample IDs must be unique. Duplicates: {duplicated}")
    return samples


def ensure_output_dirs(config: dict[str, Any]) -> dict[str, Path]:
    output = config.get("output", {})
    paths = {
        "processed_dir": Path(output.get("processed_dir", "data/processed")),
        "figures_dir": Path(output.get("figures_dir", "reports/figures")),
        "tables_dir": Path(output.get("tables_dir", "reports/tables")),
        "logs_dir": Path(output.get("logs_dir", "logs")),
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def load_sample(row: pd.Series) -> ad.AnnData:
    sample_id = str(row["sample_id"])
    sample_path = Path(row["path"])
    logging.info("Loading sample %s from %s", sample_id, sample_path)

    adata = sc.read_10x_mtx(sample_path, var_names="gene_symbols", make_unique=True)
    adata.obs["sample_id"] = sample_id
    adata.obs["condition"] = str(row["condition"])
    adata.obs["batch"] = str(row["batch"])
    adata.obs_names = [f"{sample_id}:{barcode}" for barcode in adata.obs_names]
    return adata


def load_and_concatenate(samples: pd.DataFrame) -> ad.AnnData:
    adatas = [load_sample(row) for _, row in samples.iterrows()]
    if not adatas:
        raise ValueError("Sample sheet contains no samples.")

    return ad.concat(
        adatas,
        join="outer",
        label="sample_id_from_concat",
        keys=samples["sample_id"].astype(str).tolist(),
        index_unique=None,
    )


def add_qc_annotations(adata: ad.AnnData, config: dict[str, Any]) -> None:
    qc_config = config.get("qc", {})
    mt_prefix = str(qc_config.get("mt_gene_prefix", "MT-"))
    ribo_prefixes = tuple(qc_config.get("ribo_gene_prefixes", ["RPS", "RPL"]))

    adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
    adata.var["ribo"] = adata.var_names.str.startswith(ribo_prefixes)
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )


def filtering_mask(adata: ad.AnnData, config: dict[str, Any]) -> pd.Series:
    qc_config = config.get("qc", {})
    mask = pd.Series(True, index=adata.obs_names)

    max_pct_mt = qc_config.get("max_pct_mt")
    max_genes = qc_config.get("max_genes")
    max_counts = qc_config.get("max_counts")

    if max_pct_mt is not None:
        mask &= adata.obs["pct_counts_mt"] <= float(max_pct_mt)
    if max_genes is not None:
        mask &= adata.obs["n_genes_by_counts"] <= int(max_genes)
    if max_counts is not None:
        mask &= adata.obs["total_counts"] <= float(max_counts)

    return mask


def summarize_filtering(
    before: pd.DataFrame,
    after: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    qc_config = config.get("qc", {})

    for sample_id in sorted(before["sample_id"].unique()):
        before_count = int((before["sample_id"] == sample_id).sum())
        after_count = int((after["sample_id"] == sample_id).sum())
        rows.append(
            {
                "sample_id": sample_id,
                "cells_before_filtering": before_count,
                "cells_after_filtering": after_count,
                "cells_removed": before_count - after_count,
                "min_genes": qc_config.get("min_genes"),
                "min_cells": qc_config.get("min_cells"),
                "max_pct_mt": qc_config.get("max_pct_mt"),
                "max_genes": qc_config.get("max_genes"),
                "max_counts": qc_config.get("max_counts"),
            }
        )

    rows.append(
        {
            "sample_id": "all",
            "cells_before_filtering": int(before.shape[0]),
            "cells_after_filtering": int(after.shape[0]),
            "cells_removed": int(before.shape[0] - after.shape[0]),
            "min_genes": qc_config.get("min_genes"),
            "min_cells": qc_config.get("min_cells"),
            "max_pct_mt": qc_config.get("max_pct_mt"),
            "max_genes": qc_config.get("max_genes"),
            "max_counts": qc_config.get("max_counts"),
        }
    )

    return pd.DataFrame(rows)


def apply_filtering(adata: ad.AnnData, config: dict[str, Any]) -> tuple[ad.AnnData, pd.DataFrame]:
    qc_config = config.get("qc", {})
    obs_before = adata.obs.copy()

    min_genes = qc_config.get("min_genes")
    min_cells = qc_config.get("min_cells")

    if min_genes is not None:
        sc.pp.filter_cells(adata, min_genes=int(min_genes))
    if min_cells is not None:
        sc.pp.filter_genes(adata, min_cells=int(min_cells))

    mask = filtering_mask(adata, config)
    adata = adata[mask].copy()
    summary = summarize_filtering(obs_before, adata.obs.copy(), config)
    return adata, summary


def preprocess(adata: ad.AnnData, config: dict[str, Any]) -> ad.AnnData:
    adata.layers["counts"] = adata.X.copy()

    normalization = config.get("normalization", {})
    sc.pp.normalize_total(
        adata,
        target_sum=float(normalization.get("target_sum", 10000)),
    )
    if normalization.get("log1p", True):
        sc.pp.log1p(adata)

    hvg = config.get("hvg", {})
    hvg_kwargs: dict[str, Any] = {
        "n_top_genes": int(hvg.get("n_top_genes", 2000)),
        "flavor": hvg.get("flavor", "seurat_v3"),
        "inplace": True,
    }
    if hvg.get("batch_key"):
        hvg_kwargs["batch_key"] = hvg["batch_key"]
    if hvg_kwargs["flavor"] == "seurat_v3":
        hvg_kwargs["layer"] = "counts"
    sc.pp.highly_variable_genes(adata, **hvg_kwargs)

    if "highly_variable" in adata.var:
        adata.raw = adata
        adata = adata[:, adata.var["highly_variable"]].copy()

    sc.pp.scale(adata, max_value=10)

    pca = config.get("pca", {})
    sc.tl.pca(
        adata,
        n_comps=int(pca.get("n_comps", 50)),
        svd_solver=pca.get("svd_solver", "arpack"),
    )
    return adata


def run_harmony_if_enabled(adata: ad.AnnData, config: dict[str, Any]) -> None:
    harmony = config.get("harmony", {})
    if not harmony.get("enabled", False):
        return

    batch_key = harmony.get("batch_key")
    if not batch_key:
        raise ValueError("Harmony is enabled but harmony.batch_key is missing.")
    if str(batch_key).lower() in BIOLOGICAL_HARMONY_KEYS:
        blocked = ", ".join(sorted(BIOLOGICAL_HARMONY_KEYS))
        raise ValueError(
            f"Refusing to use biological label '{batch_key}' as a Harmony batch key. "
            f"Blocked keys: {blocked}."
        )
    if batch_key not in adata.obs:
        raise ValueError(f"Harmony batch key is not present in adata.obs: {batch_key}")

    basis = harmony.get("basis", "X_pca")
    if basis not in adata.obsm:
        raise ValueError(f"Harmony basis is not present in adata.obsm: {basis}")

    unique_batches = adata.obs[batch_key].dropna().unique()
    if len(unique_batches) < 2:
        logging.warning(
            "Skipping Harmony because batch key %s has fewer than 2 unique non-null values. "
            "Neighbors will use regular X_pca.",
            batch_key,
        )
        return

    adjusted_basis = harmony.get("adjusted_basis", "X_pca_harmony")
    logging.warning(
        "Running Harmony on %s using batch key %s. Inspect results for possible overcorrection.",
        basis,
        batch_key,
    )
    sc.external.pp.harmony_integrate(
        adata,
        key=batch_key,
        basis=basis,
        adjusted_basis=adjusted_basis,
    )


def compute_embedding_and_clusters(adata: ad.AnnData, config: dict[str, Any]) -> None:
    harmony = config.get("harmony", {})
    neighbors = config.get("neighbors", {})
    use_rep = None
    if harmony.get("enabled", False):
        adjusted_basis = harmony.get("adjusted_basis", "X_pca_harmony")
        if adjusted_basis in adata.obsm:
            use_rep = adjusted_basis

    sc.pp.neighbors(
        adata,
        n_neighbors=int(neighbors.get("n_neighbors", 15)),
        n_pcs=int(neighbors.get("n_pcs", 30)),
        use_rep=use_rep,
    )

    umap = config.get("umap", {})
    sc.tl.umap(adata, random_state=int(umap.get("random_state", 0)))

    clustering = config.get("clustering", {})
    method = clustering.get("method", "leiden")
    if method != "leiden":
        raise ValueError(f"Unsupported clustering method: {method}")
    sc.tl.leiden(
        adata,
        resolution=float(clustering.get("resolution", 0.5)),
        key_added=clustering.get("key_added", "leiden"),
    )


def compute_markers(adata: ad.AnnData, config: dict[str, Any]) -> pd.DataFrame:
    markers = config.get("markers", {})
    groupby = markers.get("groupby", config.get("clustering", {}).get("key_added", "leiden"))
    if groupby not in adata.obs:
        raise ValueError(f"Marker groupby key is not present in adata.obs: {groupby}")

    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=markers.get("method", "wilcoxon"),
        n_genes=int(markers.get("n_genes", 100)),
    )
    return sc.get.rank_genes_groups_df(adata, group=None)


def save_plots(adata: ad.AnnData, paths: dict[str, Path], config: dict[str, Any]) -> None:
    figures_dir = paths["figures_dir"]
    sc.settings.figdir = str(figures_dir)

    qc_keys = ["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo"]
    sc.pl.violin(
        adata,
        keys=qc_keys,
        groupby="sample_id",
        rotation=90,
        multi_panel=True,
        show=False,
        save="_qc_by_sample.png",
    )
    sc.pl.scatter(
        adata,
        x="total_counts",
        y="pct_counts_mt",
        color="sample_id",
        show=False,
        save="_qc_counts_vs_mt.png",
    )

    clustering_key = config.get("clustering", {}).get("key_added", "leiden")
    umap_color = [clustering_key, "sample_id", "condition", "batch"]
    sc.pl.umap(
        adata,
        color=[key for key in umap_color if key in adata.obs],
        show=False,
        save="_clusters_samples_conditions.png",
    )


def output_h5ad_path(paths: dict[str, Path], config: dict[str, Any]) -> Path:
    return paths["processed_dir"] / "processed_with_harmony.h5ad"


def main() -> None:
    setup_logging()
    args = parse_args()

    config = read_yaml(args.config)
    sample_sheet = resolve_sample_sheet(args.config, config, args.sample_sheet)
    samples = read_samples(sample_sheet)
    paths = ensure_output_dirs(config)

    adata = load_and_concatenate(samples)
    add_qc_annotations(adata, config)
    adata, qc_summary = apply_filtering(adata, config)
    adata = preprocess(adata, config)
    run_harmony_if_enabled(adata, config)
    compute_embedding_and_clusters(adata, config)
    marker_table = compute_markers(adata, config)

    adata.uns["analysis_params"] = {
        "config": config,
        "config_path": str(args.config),
        "sample_sheet": str(sample_sheet),
    }

    qc_summary.to_csv(paths["tables_dir"] / "qc_filtering_summary.csv", index=False)
    marker_table.to_csv(paths["tables_dir"] / "cluster_markers.csv", index=False)
    save_plots(adata, paths, config)
    adata.write_h5ad(output_h5ad_path(paths, config))


if __name__ == "__main__":
    main()

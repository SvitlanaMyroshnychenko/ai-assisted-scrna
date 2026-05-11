#!/usr/bin/env python
"""Config-driven scRNA-seq workflow for remote analysis.

This script is intended to be executed on the remote analysis server where the
real sample sheet and input count matrices are available. The local repository
contains templates only.
"""

from __future__ import annotations

import argparse
from datetime import datetime
import logging
from pathlib import Path
from typing import Any

import anndata as ad
import harmonypy as hp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import yaml


DEFAULT_CONFIG = Path("configs/pipeline.template.yaml")
DEFAULT_TEMPLATE_SAMPLE_SHEET = Path("configs/samples.template.csv")
REQUIRED_SAMPLE_COLUMNS = ("sample_id", "path", "condition", "batch")
RISKY_HARMONY_BATCH_KEY_TERMS = (
    "disease",
    "diagnosis",
    "clinical",
    "pathology",
    "pathological",
    "braak",
    "cerad",
    "cog",
    "mmse",
    "phenotype",
    "treatment",
    "outcome",
    "sex",
    "age",
)
RISKY_HARMONY_BATCH_KEY_WARNING = (
    "This batch key may encode biological or clinical variation and should not be used "
    "for Harmony unless explicitly justified."
)
STEPS = ("qc", "preprocess", "harmony", "cluster", "markers", "plots")
STEP_H5AD_FILES = {
    "qc": "qc.h5ad",
    "preprocess": "preprocessed.h5ad",
    "harmony": "harmony.h5ad",
    "cluster": "clustered.h5ad",
    "markers": "markers.h5ad",
    "plots": "plots.h5ad",
}


def input_format(config: dict[str, Any]) -> str:
    return str(config.get("input", {}).get("format", "10x_mtx")).lower()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a reproducible scRNA-seq analysis workflow."
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
    parser.add_argument(
        "--step",
        choices=STEPS,
        default=None,
        help="Run one pipeline step and save an intermediate h5ad under results/<run_id>/data/.",
    )
    parser.add_argument(
        "--plot-context",
        choices=("auto", "qc", "preprocess", "harmony", "cluster", "markers"),
        default="auto",
        help="Select which step intermediate to use for --step plots.",
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


def resolve_run_id(config: dict[str, Any], step: str | None = None) -> str:
    run_config = config.get("run", {}) or {}
    run_id = str(run_config.get("id", config.get("project_name", "scrna_run")))
    if not run_config.get("add_timestamp", False):
        return run_id

    base_dir = Path(run_config.get("base_dir", "results"))
    if step and step != "qc":
        existing_runs = sorted(
            [path for path in base_dir.glob(f"{run_id}_*") if path.is_dir()],
            key=lambda path: path.name,
        )
        if not existing_runs:
            raise FileNotFoundError(
                f"No timestamped run found for run.id '{run_id}' in {base_dir}. "
                "Run the 'qc' step first or set run.add_timestamp=false with a fixed run.id."
            )
        return existing_runs[-1].name

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    candidate = f"{run_id}_{timestamp}"
    suffix = 1
    while (base_dir / candidate).exists():
        candidate = f"{run_id}_{timestamp}_{suffix:02d}"
        suffix += 1
    return candidate


def ensure_run_output_dirs(config: dict[str, Any], run_id: str) -> dict[str, Path]:
    base_dir = Path((config.get("run", {}) or {}).get("base_dir", "results")) / run_id
    paths = {
        "run_dir": base_dir,
        "processed_dir": base_dir / "data",
        "figures_dir": base_dir / "figures",
        "tables_dir": base_dir / "tables",
        "logs_dir": base_dir / "logs",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def intermediate_h5ad_path(paths: dict[str, Path], step: str) -> Path:
    return paths["processed_dir"] / STEP_H5AD_FILES[step]


def read_intermediate(paths: dict[str, Path], step: str) -> ad.AnnData:
    path = intermediate_h5ad_path(paths, step)
    if not path.exists():
        raise FileNotFoundError(
            f"Missing intermediate for step '{step}': {path}. "
            f"Run the '{step}' step first."
        )
    logging.info("Loading intermediate %s", path)
    return sc.read_h5ad(path)


def write_intermediate(adata: ad.AnnData, paths: dict[str, Path], step: str) -> None:
    path = intermediate_h5ad_path(paths, step)
    logging.info("Writing %s intermediate to %s", step, path)
    adata.write_h5ad(path)


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


def load_h5ad(config: dict[str, Any]) -> ad.AnnData:
    input_config = config.get("input", {})
    input_path = input_config.get("path")
    if not input_path:
        raise ValueError("Missing input.path in config for input.format: h5ad.")

    path = Path(input_path)
    logging.info("Loading h5ad input from %s", path)
    return sc.read_h5ad(path)


def load_input(
    config_path: Path,
    config: dict[str, Any],
    sample_sheet_arg: Path | None,
) -> tuple[ad.AnnData, Path | None]:
    fmt = input_format(config)
    if fmt == "10x_mtx":
        sample_sheet = resolve_sample_sheet(config_path, config, sample_sheet_arg)
        samples = read_samples(sample_sheet)
        return load_and_concatenate(samples), sample_sheet
    if fmt == "h5ad":
        if sample_sheet_arg is not None:
            logging.warning("Ignoring --sample-sheet because input.format is h5ad.")
        return load_h5ad(config), None
    raise ValueError(f"Unsupported input.format: {fmt}")


def apply_debug_subset(adata: ad.AnnData, config: dict[str, Any]) -> ad.AnnData:
    debug = config.get("debug", {}) or {}
    max_cells = debug.get("max_cells")
    if max_cells is None:
        return adata

    max_cells = int(max_cells)
    if max_cells <= 0:
        raise ValueError("debug.max_cells must be a positive integer when set.")

    selected_count = min(max_cells, adata.n_obs)
    random_state = int(debug.get("random_state", 42))
    selected_obs_names = adata.obs.sample(
        n=selected_count,
        random_state=random_state,
    ).index
    logging.info(
        "Debug mode selected %s of %s cells using random_state=%s.",
        selected_count,
        adata.n_obs,
        random_state,
    )
    return adata[selected_obs_names].copy()


def add_qc_annotations(adata: ad.AnnData, config: dict[str, Any]) -> None:
    qc_config = config.get("qc", {})
    mt_prefix = str(qc_config.get("mt_gene_prefix", "MT-"))
    ribo_prefixes = tuple(qc_config.get("ribo_gene_prefixes", ["RPS", "RPL"]))
    hb_gene_pattern = qc_config.get("hb_gene_pattern", r"^HB[^(P)]")

    adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
    adata.var["ribo"] = adata.var_names.str.startswith(ribo_prefixes)
    qc_vars = ["mt", "ribo"]
    if hb_gene_pattern:
        adata.var["hb"] = adata.var_names.str.contains(str(hb_gene_pattern), regex=True)
        qc_vars.append("hb")
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=qc_vars,
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    adata.obs["log1p_total_counts"] = np.log1p(adata.obs["total_counts"])
    adata.obs["log1p_n_genes_by_counts"] = np.log1p(adata.obs["n_genes_by_counts"])


def mad_outlier_mask(values: pd.Series, nmads: float, direction: str = "two_sided") -> pd.Series:
    median = values.median()
    mad = (values - median).abs().median()
    if pd.isna(mad) or mad == 0:
        return pd.Series(False, index=values.index)

    lower = median - (nmads * mad)
    upper = median + (nmads * mad)
    if direction == "upper":
        return values > upper
    if direction != "two_sided":
        raise ValueError(f"Unsupported MAD outlier direction: {direction}")
    return (values < lower) | (values > upper)


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

    mad_filter = qc_config.get("mad_filter", {}) or {}
    if mad_filter.get("enabled", False):
        mt_direction = str(mad_filter.get("mt_direction", "upper"))
        if mt_direction not in {"upper", "two_sided"}:
            raise ValueError(
                "qc.mad_filter.mt_direction must be one of: upper, two_sided. "
                f"Found: {mt_direction}"
            )
        mad_fields = {
            "log1p_total_counts": mad_filter.get("log1p_total_counts_nmads", 5),
            "log1p_n_genes_by_counts": mad_filter.get("log1p_n_genes_by_counts_nmads", 5),
        }
        for column, nmads in mad_fields.items():
            if nmads is None:
                continue
            if column not in adata.obs:
                raise ValueError(f"MAD QC column is not present in adata.obs: {column}")
            mask &= ~mad_outlier_mask(adata.obs[column], float(nmads))
        mt_nmads = mad_filter.get("pct_counts_mt_nmads", 3)
        if mt_nmads is not None:
            if "pct_counts_mt" not in adata.obs:
                raise ValueError("MAD QC column is not present in adata.obs: pct_counts_mt")
            mask &= ~mad_outlier_mask(adata.obs["pct_counts_mt"], float(mt_nmads), mt_direction)

    return mask


def summarize_filtering(
    before: pd.DataFrame,
    after: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    qc_config = config.get("qc", {})
    mad_filter = qc_config.get("mad_filter", {}) or {}
    qc_params = {
        "min_genes": qc_config.get("min_genes"),
        "min_cells": qc_config.get("min_cells"),
        "max_pct_mt": qc_config.get("max_pct_mt"),
        "max_genes": qc_config.get("max_genes"),
        "max_counts": qc_config.get("max_counts"),
        "mad_filter_enabled": mad_filter.get("enabled", False),
        "log1p_total_counts_nmads": mad_filter.get("log1p_total_counts_nmads", 5),
        "log1p_n_genes_by_counts_nmads": mad_filter.get(
            "log1p_n_genes_by_counts_nmads", 5
        ),
        "pct_counts_mt_nmads": mad_filter.get("pct_counts_mt_nmads", 3),
        "mt_direction": mad_filter.get("mt_direction", "upper"),
    }

    if "sample_id" in before:
        for sample_id in sorted(before["sample_id"].dropna().unique()):
            before_count = int((before["sample_id"] == sample_id).sum())
            after_count = int((after["sample_id"] == sample_id).sum())
            rows.append(
                {
                    "sample_id": sample_id,
                    "cells_before_filtering": before_count,
                    "cells_after_filtering": after_count,
                    "cells_removed": before_count - after_count,
                    **qc_params,
                }
            )

    rows.append(
        {
            "sample_id": "all",
            "cells_before_filtering": int(before.shape[0]),
            "cells_after_filtering": int(after.shape[0]),
            "cells_removed": int(before.shape[0] - after.shape[0]),
            **qc_params,
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
    adata.raw = adata.copy()

    hvg = config.get("hvg", {})
    hvg_kwargs: dict[str, Any] = {
        "n_top_genes": int(hvg.get("n_top_genes", 3000)),
        "flavor": hvg.get("flavor", "seurat_v3"),
        "inplace": True,
        "subset": False,
    }
    if hvg.get("batch_key"):
        if hvg["batch_key"] not in adata.obs:
            raise ValueError(f"HVG batch key is not present in adata.obs: {hvg['batch_key']}")
        hvg_kwargs["batch_key"] = hvg["batch_key"]
    if hvg_kwargs["flavor"] == "seurat_v3":
        hvg_kwargs["layer"] = "counts"
    sc.pp.highly_variable_genes(adata, **hvg_kwargs)

    scaling = config.get("scaling", {}) or {}
    sc.pp.scale(adata, max_value=float(scaling.get("max_value", 10)))

    pca = config.get("pca", {})
    sc.tl.pca(
        adata,
        n_comps=int(pca.get("n_comps", 50)),
        use_highly_variable=True,
        svd_solver=pca.get("svd_solver", "arpack"),
        random_state=int(pca.get("random_state", 42)),
    )
    return adata


def _shape_text(values: Any | None) -> str:
    return str(tuple(values.shape)) if values is not None and hasattr(values, "shape") else ""


def _batch_groups_preview(values: pd.Series, max_groups: int = 10) -> str:
    groups = pd.Series(values.dropna().unique()).astype(str).tolist()
    preview = groups[:max_groups]
    suffix = f"; ... ({len(groups)} total)" if len(groups) > max_groups else ""
    return "; ".join(preview) + suffix


def _risky_harmony_batch_key_warning(batch_key: str | None) -> str:
    if not batch_key:
        return ""
    normalized = str(batch_key).lower()
    if any(term in normalized for term in RISKY_HARMONY_BATCH_KEY_TERMS):
        return RISKY_HARMONY_BATCH_KEY_WARNING
    return ""


def _harmony_summary_row(
    *,
    harmony_enabled: bool,
    harmony_ran: bool,
    skip_reason: str,
    batch_key: str,
    n_batch_groups: int | None,
    batch_groups_preview: str,
    basis: str,
    basis_shape: str,
    adjusted_basis: str,
    adjusted_basis_shape: str,
    risky_batch_key_warning: str,
) -> dict[str, Any]:
    return {
        "harmony_enabled": harmony_enabled,
        "harmony_ran": harmony_ran,
        "harmony_skipped": not harmony_ran,
        "skip_reason": skip_reason,
        "batch_key": batch_key,
        "n_batch_groups": n_batch_groups,
        "batch_groups_preview": batch_groups_preview,
        "basis": basis,
        "basis_shape": basis_shape,
        "adjusted_basis": adjusted_basis,
        "adjusted_basis_shape": adjusted_basis_shape,
        "risky_batch_key_warning": risky_batch_key_warning,
    }


def run_harmony_if_enabled(adata: ad.AnnData, config: dict[str, Any]) -> pd.DataFrame:
    harmony = config.get("harmony", {}) or {}
    harmony_enabled = bool(harmony.get("enabled", False))
    batch_key = str(harmony.get("batch_key") or "")
    basis = str(harmony.get("basis") or "X_pca")
    adjusted_basis = str(harmony.get("adjusted_basis") or "X_pca_harmony")
    basis_values = adata.obsm.get(basis)
    risky_warning = _risky_harmony_batch_key_warning(batch_key)

    if risky_warning:
        logging.warning(risky_warning)
        print(f"WARNING: {risky_warning}")

    if not harmony_enabled:
        return pd.DataFrame(
            [
                _harmony_summary_row(
                    harmony_enabled=False,
                    harmony_ran=False,
                    skip_reason="Harmony skipped because harmony.enabled is false.",
                    batch_key=batch_key,
                    n_batch_groups=None,
                    batch_groups_preview="",
                    basis=basis,
                    basis_shape=_shape_text(basis_values),
                    adjusted_basis=adjusted_basis,
                    adjusted_basis_shape=_shape_text(adata.obsm.get(adjusted_basis)),
                    risky_batch_key_warning=risky_warning,
                )
            ]
        )

    if not batch_key:
        raise ValueError("Harmony is enabled but harmony.batch_key is missing or empty.")
    if batch_key not in adata.obs:
        available_columns = list(map(str, adata.obs.columns[:25]))
        raise ValueError(
            f"Harmony batch key '{batch_key}' is not present in adata.obs. "
            f"Available obs columns preview: {available_columns}"
        )

    batch_values = adata.obs[batch_key]
    unique_batches = batch_values.dropna().unique()
    n_batch_groups = len(unique_batches)
    batch_preview = _batch_groups_preview(batch_values)
    if n_batch_groups < 2:
        skip_reason = (
            f"Harmony skipped because batch key '{batch_key}' has fewer than 2 "
            "non-null groups."
        )
        logging.warning(skip_reason)
        return pd.DataFrame(
            [
                _harmony_summary_row(
                    harmony_enabled=True,
                    harmony_ran=False,
                    skip_reason=skip_reason,
                    batch_key=batch_key,
                    n_batch_groups=n_batch_groups,
                    batch_groups_preview=batch_preview,
                    basis=basis,
                    basis_shape=_shape_text(basis_values),
                    adjusted_basis=adjusted_basis,
                    adjusted_basis_shape=_shape_text(adata.obsm.get(adjusted_basis)),
                    risky_batch_key_warning=risky_warning,
                )
            ]
        )

    if basis not in adata.obsm:
        raise ValueError(f"Harmony basis '{basis}' is missing from adata.obsm.")

    if basis_values.ndim != 2:
        raise ValueError(
            f"Harmony basis {basis} must be 2-dimensional. Found ndim={basis_values.ndim}."
        )
    if basis_values.shape[0] != adata.n_obs:
        raise ValueError(
            f"Harmony basis {basis} has shape {basis_values.shape}; expected first "
            f"dimension to match adata.n_obs={adata.n_obs}."
        )

    logging.warning(
        "Running Harmony on %s using batch key %s. Inspect results for possible overcorrection.",
        basis,
        batch_key,
    )
    harmony_out = hp.run_harmony(basis_values, adata.obs, batch_key)
    corrected = np.asarray(harmony_out.Z_corr)
    if corrected.shape == basis_values.shape:
        adjusted_values = corrected
    elif corrected.T.shape == basis_values.shape:
        adjusted_values = corrected.T
    else:
        raise ValueError(
            f"Harmony output shape mismatch. PCA basis '{basis}' has shape "
            f"{basis_values.shape}; Harmony returned shape {corrected.shape}."
        )

    adata.obsm[adjusted_basis] = adjusted_values

    return pd.DataFrame(
        [
            _harmony_summary_row(
                harmony_enabled=True,
                harmony_ran=True,
                skip_reason="",
                batch_key=batch_key,
                n_batch_groups=n_batch_groups,
                batch_groups_preview=batch_preview,
                basis=basis,
                basis_shape=_shape_text(basis_values),
                adjusted_basis=adjusted_basis,
                adjusted_basis_shape=_shape_text(adata.obsm.get(adjusted_basis)),
                risky_batch_key_warning=risky_warning,
            )
        ]
    )


def choose_neighbors_rep(adata: ad.AnnData, config: dict[str, Any]) -> str:
    harmony = config.get("harmony", {}) or {}
    adjusted_basis = str(harmony.get("adjusted_basis") or "X_pca_harmony")
    if harmony.get("enabled", False) and adjusted_basis in adata.obsm:
        return adjusted_basis
    return "X_pca"


def _cluster_sizes_preview(adata: ad.AnnData, clustering_key: str, max_clusters: int = 10) -> str:
    if clustering_key not in adata.obs:
        return ""
    sizes = adata.obs[clustering_key].astype(str).value_counts().sort_index()
    preview = [f"{cluster}:{count}" for cluster, count in sizes.head(max_clusters).items()]
    suffix = f"; ... ({len(sizes)} total)" if len(sizes) > max_clusters else ""
    return "; ".join(preview) + suffix


def compute_embedding_and_clusters(
    adata: ad.AnnData,
    config: dict[str, Any],
    input_h5ad: Path | None = None,
) -> pd.DataFrame:
    neighbors = config.get("neighbors", {}) or {}
    use_rep = choose_neighbors_rep(adata, config)
    if use_rep not in adata.obsm:
        raise ValueError(f"Clustering basis '{use_rep}' is missing from adata.obsm.")

    sc.pp.neighbors(
        adata,
        n_neighbors=int(neighbors.get("n_neighbors", 15)),
        n_pcs=int(neighbors.get("n_pcs", 30)),
        use_rep=use_rep,
    )

    umap = config.get("umap", {})
    sc.tl.umap(adata, random_state=int(umap.get("random_state", 42)))

    clustering = config.get("clustering", {})
    method = clustering.get("method", "leiden")
    if method != "leiden":
        raise ValueError(f"Unsupported clustering method: {method}")
    sc.tl.leiden(
        adata,
        resolution=float(clustering.get("resolution", 0.5)),
        key_added=clustering.get("key_added", "leiden"),
        random_state=int(clustering.get("random_state", 42)),
    )

    clustering_key = str(clustering.get("key_added", "leiden"))
    x_umap = adata.obsm.get("X_umap")
    cluster_sizes = _cluster_sizes_preview(adata, clustering_key)
    n_clusters = (
        int(adata.obs[clustering_key].nunique(dropna=True))
        if clustering_key in adata.obs
        else 0
    )
    return pd.DataFrame(
        [
            {
                "input_h5ad": str(input_h5ad) if input_h5ad is not None else "",
                "cells": adata.n_obs,
                "genes": adata.n_vars,
                "neighbors_use_rep": use_rep,
                "n_neighbors": int(neighbors.get("n_neighbors", 15)),
                "n_pcs": int(neighbors.get("n_pcs", 30)),
                "umap_random_state": int(umap.get("random_state", 42)),
                "clustering_method": method,
                "leiden_resolution": float(clustering.get("resolution", 0.5)),
                "clustering_key": clustering_key,
                "n_clusters": n_clusters,
                "cluster_sizes_preview": cluster_sizes,
                "X_umap_present": x_umap is not None,
                "X_umap_shape": _shape_text(x_umap),
                "X_pca_present": "X_pca" in adata.obsm,
                "X_pca_harmony_present": "X_pca_harmony" in adata.obsm,
            }
        ]
    )


def compute_markers(
    adata: ad.AnnData,
    config: dict[str, Any],
    input_h5ad: Path | None = None,
    marker_table_path: Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    markers = config.get("markers", {}) or {}
    groupby = str(markers.get("groupby", config.get("clustering", {}).get("key_added", "leiden")))
    method = str(markers.get("method", "wilcoxon"))
    n_genes = int(markers.get("n_genes", 100))
    if groupby not in adata.obs:
        raise ValueError(f"Marker groupby key is not present in adata.obs: {groupby}")

    groups = pd.Series(adata.obs[groupby].dropna().unique()).astype(str).tolist()
    groups = sorted(groups)
    if len(groups) < 2:
        raise ValueError(
            f"Marker groupby key '{groupby}' must contain at least 2 non-null groups. "
            f"Found {len(groups)}."
        )

    use_raw = adata.raw is not None
    if not use_raw:
        logging.warning("adata.raw is missing; marker testing will use current adata.X.")
        print("WARNING: adata.raw is missing; marker testing will use current adata.X.")

    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method,
        n_genes=n_genes,
        use_raw=use_raw,
    )

    marker_table = sc.get.rank_genes_groups_df(adata, group=None).rename(
        columns={
            "group": "cluster",
            "names": "gene",
            "logfoldchanges": "logfoldchange",
            "pvals": "pval",
            "pvals_adj": "pval_adj",
        }
    )
    marker_table["rank"] = marker_table.groupby("cluster").cumcount() + 1
    marker_table = marker_table[
        ["cluster", "gene", "rank", "scores", "logfoldchange", "pval", "pval_adj"]
    ].rename(columns={"scores": "score"})

    markers_summary = pd.DataFrame(
        [
            {
                "input_h5ad": str(input_h5ad) if input_h5ad is not None else "",
                "cells": adata.n_obs,
                "genes": adata.n_vars,
                "groupby": groupby,
                "method": method,
                "n_groups": len(groups),
                "groups": "; ".join(groups),
                "n_genes_requested": n_genes,
                "used_raw": use_raw,
                "marker_table_path": str(marker_table_path) if marker_table_path is not None else "",
            }
        ]
    )
    return marker_table, markers_summary


def configure_plot_style() -> None:
    sc.set_figure_params(dpi=120, color_map="viridis", frameon=False)
    try:
        import seaborn as sns

        sns.set_style("whitegrid")
    except ImportError:
        logging.warning("seaborn is not installed; continuing without seaborn whitegrid style.")


def configure_qc_plot_style() -> None:
    plt.rcParams.update(plt.rcParamsDefault)
    sc.set_figure_params(dpi=100, color_map="viridis", frameon=False)
    try:
        import seaborn as sns

        sns.set_style("whitegrid")
    except ImportError:
        logging.warning("seaborn is not installed; continuing without seaborn whitegrid style.")


def save_preprocess_plots(adata: ad.AnnData, figures_dir: Path) -> None:
    if "highly_variable" in adata.var:
        try:
            configure_qc_plot_style()
            with plt.rc_context({"figure.figsize": (8, 4)}):
                sc.pl.highly_variable_genes(adata, show=False)
            fig = plt.gcf()
            fig.suptitle("Highly variable genes", y=1.02)
            fig.tight_layout()
            fig.savefig(figures_dir / "hvg_plot.png", bbox_inches="tight")
            plt.close(fig)
        except Exception as exc:
            plt.close("all")
            print(f"WARNING: skipping hvg_plot.png because HVG plotting failed: {exc}")
    else:
        print("WARNING: skipping hvg_plot.png because highly_variable is not present in adata.var.")

    if "pca" in adata.uns and "variance_ratio" in adata.uns["pca"]:
        configure_qc_plot_style()
        variance_ratio = np.asarray(adata.uns["pca"]["variance_ratio"])
        pcs = np.arange(1, len(variance_ratio) + 1)
        with plt.rc_context({"figure.figsize": (8, 5), "font.size": 10}):
            fig, ax1 = plt.subplots(dpi=100)
            ax1.bar(
                pcs,
                variance_ratio,
                color="steelblue",
                alpha=0.8,
                label="Individual Variance",
            )
            ax1.axvline(
                30,
                color="black",
                linestyle=":",
                linewidth=1.5,
                label="PC 30 Threshold",
            )
            ax1.set_title("PCA Scree Plot", fontsize=14)
            ax1.set_xlabel("Principal Components (PCs)")
            ax1.set_ylabel("Explained Variance Ratio")
            ax1.grid(linestyle="--", alpha=0.5)
            ax2 = ax1.twinx()
            ax2.plot(
                pcs,
                np.cumsum(variance_ratio),
                color="firebrick",
                marker="o",
                markersize=4,
                label="Cumulative Variance",
            )
            ax2.set_ylabel("Cumulative Variance Ratio")
            lines1, labels1 = ax1.get_legend_handles_labels()
            lines2, labels2 = ax2.get_legend_handles_labels()
            ax1.legend(lines1 + lines2, labels1 + labels2, loc="center right", frameon=False)
            fig.tight_layout()
            fig.savefig(figures_dir / "pca_scree.png", bbox_inches="tight")
            plt.close(fig)
    else:
        print("WARNING: skipping pca_scree.png because PCA variance_ratio is not present.")

    required_obs = ["total_counts", "pct_counts_mt"]
    missing_obs = [key for key in required_obs if key not in adata.obs]
    if "X_pca" not in adata.obsm:
        print("WARNING: skipping pca_scatter_qc_metrics.png because X_pca is not present in adata.obsm.")
        return
    if adata.obsm["X_pca"].shape[1] < 2:
        print("WARNING: skipping pca_scatter_qc_metrics.png because X_pca has fewer than 2 components.")
        return
    if missing_obs:
        print(
            "WARNING: skipping pca_scatter_qc_metrics.png because required obs "
            f"columns are missing: {missing_obs}"
        )
        return

    configure_qc_plot_style()
    pca_values = adata.obsm["X_pca"]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), dpi=100)
    for ax, color_key, title in zip(
        axes,
        required_obs,
        ["PC1 vs PC2 colored by total counts", "PC1 vs PC2 colored by mitochondrial percentage"],
    ):
        scatter = ax.scatter(
            pca_values[:, 0],
            pca_values[:, 1],
            c=adata.obs[color_key],
            s=10,
            alpha=0.6,
            cmap="viridis",
            edgecolors="none",
            linewidths=0,
        )
        ax.set_title(title)
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        for spine in ax.spines.values():
            spine.set_visible(False)
        fig.colorbar(scatter, ax=ax, label=color_key, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(figures_dir / "pca_scatter_qc_metrics.png", bbox_inches="tight")
    plt.close(fig)


def has_preprocess_outputs(adata: ad.AnnData) -> bool:
    return (
        "highly_variable" in adata.var
        or "pca" in adata.uns
        or "X_pca" in adata.obsm
    )


def save_scanpy_umap(
    adata: ad.AnnData,
    figures_dir: Path,
    filename: str,
    *,
    color: str | list[str],
    title: str | list[str],
    palette: Any | None = None,
    ncols: int | None = None,
    size: float = 3,
    legend_loc: str | None = "right margin",
) -> None:
    kwargs: dict[str, Any] = {
        "color": color,
        "title": title,
        "size": size,
        "frameon": False,
        "show": False,
    }
    if palette is not None:
        kwargs["palette"] = palette
    if ncols is not None:
        kwargs["ncols"] = ncols
    if legend_loc is not None:
        kwargs["legend_loc"] = legend_loc

    sc.pl.umap(adata, **kwargs)
    fig = plt.gcf()
    fig.tight_layout()
    fig.savefig(figures_dir / filename, bbox_inches="tight")
    plt.close(fig)


def save_clustering_plots(adata: ad.AnnData, paths: dict[str, Path], config: dict[str, Any]) -> None:
    figures_dir = paths["figures_dir"]
    configure_qc_plot_style()

    if "X_umap" not in adata.obsm:
        logging.warning("Skipping clustering UMAP plots because X_umap is not present in adata.obsm.")
        return

    clustering_key = "leiden"
    leiden_palette = sc.pl.palettes.godsnot_102

    qc_counts_genes = ["total_counts", "n_genes_by_counts"]
    missing_qc = [column for column in qc_counts_genes if column not in adata.obs]
    if missing_qc:
        logging.warning(
            "Skipping umap_qc_counts_genes.png because required obs columns are missing: %s",
            missing_qc,
        )
        print(
            "WARNING: skipping umap_qc_counts_genes.png because required obs "
            f"columns are missing: {missing_qc}"
        )
    else:
        with plt.rc_context({"figure.figsize": (10, 6)}):
            sc.pl.umap(
                adata,
                color=qc_counts_genes,
                cmap="magma",
                size=3,
                ncols=2,
                frameon=False,
                title=["Total Counts", "Number of Genes"],
                show=False,
            )
        fig = plt.gcf()
        fig.savefig(figures_dir / "umap_qc_counts_genes.png", bbox_inches="tight")
        plt.close(fig)

    if clustering_key not in adata.obs:
        logging.warning("Skipping Leiden UMAP plots because clustering key %s is missing.", clustering_key)
        print(f"WARNING: skipping Leiden UMAP plots because clustering key is missing: {clustering_key}")
    else:
        adata.uns[f"{clustering_key}_colors"] = leiden_palette
        with plt.rc_context(
            {
                "figure.figsize": (10, 8),
                "font.size": 10,
                "legend.fontsize": 10,
            }
        ):
            sc.pl.umap(
                adata,
                color=clustering_key,
                size=2,
                frameon=False,
                title="Leiden Clusters",
                legend_loc="right margin",
                legend_fontoutline=2,
                legend_fontweight="bold",
                add_outline=True,
                na_in_legend=False,
                show=False,
            )
        fig = plt.gcf()
        fig.savefig(figures_dir / "umap_leiden.png", bbox_inches="tight")
        plt.close(fig)

        sc.pl.umap(
            adata,
            color=clustering_key,
            legend_loc="on data",
            legend_fontsize=10,
            legend_fontoutline=2,
            frameon=False,
            add_outline=True,
            title="Leiden Clusters (Numbered)",
            show=False,
        )
        fig = plt.gcf()
        fig.savefig(figures_dir / "umap_leiden_numbered.png", bbox_inches="tight")
        plt.close(fig)

    dataset_label = str(config.get("dataset_id") or config.get("project_name") or "Dataset")
    target_metadata = ["leiden", "major.celltype", "Pathologic_diagnosis_of_AD", "braaksc"]
    with plt.rc_context({"figure.dpi": 120, "font.size": 14}):
        fig, axes = plt.subplots(2, 2, figsize=(16, 10))
        axes = axes.flatten()
        plotted_any = False
        for i, metadata in enumerate(target_metadata):
            if metadata not in adata.obs.columns:
                logging.warning(
                    "Skipping metadata panel for %s because it is missing from adata.obs.",
                    metadata,
                )
                print(
                    "WARNING: skipping metadata panel for "
                    f"{metadata} because it is missing from adata.obs."
                )
                axes[i].set_visible(False)
                continue

            is_leiden = metadata == "leiden"
            cmap = "YlOrRd" if metadata == "braaksc" else "magma"
            sc.pl.umap(
                adata,
                color=metadata,
                ax=axes[i],
                show=False,
                frameon=False,
                title=f"{dataset_label} dataset: {metadata}",
                add_outline=True,
                outline_width=(0.1, 0.05),
                size=3,
                palette=leiden_palette if is_leiden else None,
                cmap=cmap,
                legend_loc="on data" if is_leiden else "right margin",
                legend_fontsize=10,
                legend_fontoutline=2,
            )
            plotted_any = True

        if plotted_any:
            plt.tight_layout(pad=1.0)
            plt.subplots_adjust(wspace=0.1, hspace=0.15)
            fig.savefig(figures_dir / "umap_metadata_panel.png", bbox_inches="tight")
        else:
            logging.warning("Skipping umap_metadata_panel.png because none of the target metadata columns exist.")
            print(
                "WARNING: skipping umap_metadata_panel.png because none of the "
                "target metadata columns exist."
            )
        plt.close(fig)


def save_marker_plots(adata: ad.AnnData, paths: dict[str, Path], config: dict[str, Any]) -> None:
    figures_dir = paths["figures_dir"]
    configure_qc_plot_style()

    markers = config.get("markers", {}) or {}
    groupby = str(markers.get("groupby", config.get("clustering", {}).get("key_added", "leiden")))
    if groupby not in adata.obs:
        logging.warning("Skipping marker plots because groupby key %s is missing.", groupby)
        print(f"WARNING: skipping marker plots because groupby key is missing: {groupby}")
        return
    if "rank_genes_groups" not in adata.uns:
        logging.warning("Skipping marker plots because adata.uns['rank_genes_groups'] is missing.")
        print("WARNING: skipping marker plots because rank_genes_groups results are missing.")
        return

    try:
        with plt.rc_context({"figure.figsize": (12, 8), "font.size": 10}):
            sc.pl.rank_genes_groups_dotplot(
                adata,
                n_genes=5,
                groupby=groupby,
                show=False,
            )
        fig = plt.gcf()
        fig.savefig(figures_dir / "marker_dotplot_top_genes.png", bbox_inches="tight")
        plt.close(fig)
    except Exception as exc:
        plt.close("all")
        logging.warning("Skipping marker_dotplot_top_genes.png because plotting failed: %s", exc)
        print(f"WARNING: skipping marker_dotplot_top_genes.png because plotting failed: {exc}")


def save_plots(
    adata: ad.AnnData,
    paths: dict[str, Path],
    config: dict[str, Any],
    plot_context: str = "auto",
) -> None:
    figures_dir = paths["figures_dir"]
    configure_qc_plot_style()
    sc.settings.figdir = str(figures_dir)

    if plot_context == "cluster":
        save_clustering_plots(adata, paths, config)
        return
    if plot_context == "markers":
        save_marker_plots(adata, paths, config)
        return

    effective_context = "auto" if plot_context == "harmony" else plot_context
    should_plot_qc = effective_context in {"auto", "qc"} and "total_counts" in adata.obs
    should_plot_preprocess = effective_context in {"auto", "preprocess"} and has_preprocess_outputs(adata)
    should_plot_umap = effective_context == "auto"

    if should_plot_qc:
        qc_keys = ["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo"]
        if "pct_counts_hb" in adata.obs:
            qc_keys.append("pct_counts_hb")
        if "sample_id" in adata.obs:
            sc.pl.violin(
                adata,
                keys=qc_keys,
                groupby="sample_id",
                rotation=90,
                multi_panel=True,
                show=False,
                save="_qc_by_sample.png",
            )
        else:
            sc.pl.violin(
                adata,
                keys=qc_keys,
                rotation=90,
                multi_panel=True,
                show=False,
                save="_qc.png",
            )

        scatter_kwargs: dict[str, Any] = {
            "x": "total_counts",
            "y": "pct_counts_mt",
            "title": "Library size vs mitochondrial percentage",
            "show": False,
            "save": "_qc_counts_vs_mt.png",
        }
        if "sample_id" in adata.obs:
            scatter_kwargs["color"] = "sample_id"
        sc.pl.scatter(adata, **scatter_kwargs)

        qc_scatter_keys = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
        missing_qc_scatter_keys = [key for key in qc_scatter_keys if key not in adata.obs]
        if missing_qc_scatter_keys:
            print(
                "WARNING: skipping qc_scatter_counts_vs_genes_mt.png because required obs "
                f"columns are missing: {missing_qc_scatter_keys}"
            )
        else:
            fig, ax = plt.subplots(figsize=(5, 4), dpi=100)
            scatter = ax.scatter(
                adata.obs["total_counts"],
                adata.obs["n_genes_by_counts"],
                c=adata.obs["pct_counts_mt"],
                s=20,
                alpha=0.6,
                cmap="viridis",
                edgecolors="none",
            )
            ax.set_title("Gene diversity colored by mitochondrial percentage")
            ax.set_xlabel("Total counts")
            ax.set_ylabel("Detected genes")
            for spine in ax.spines.values():
                spine.set_visible(False)
            fig.colorbar(scatter, ax=ax, label="pct_counts_mt")
            fig.tight_layout()
            fig.savefig(figures_dir / "qc_scatter_counts_vs_genes_mt.png")
            plt.close(fig)

        qc_titles = {
            "total_counts": "Total counts per cell",
            "n_genes_by_counts": "Detected genes per cell",
            "pct_counts_mt": "Mitochondrial percentage",
            "pct_counts_ribo": "Ribosomal percentage",
            "pct_counts_hb": "Hemoglobin percentage",
        }
        for key in qc_keys:
            if key not in adata.obs:
                continue
            fig, ax = plt.subplots(figsize=(5, 4), dpi=100)
            adata.obs[key].hist(ax=ax, bins=50)
            ax.set_title(qc_titles.get(key, key))
            ax.set_xlabel(qc_titles.get(key, key))
            ax.set_ylabel("Cells")
            fig.tight_layout()
            fig.savefig(figures_dir / f"qc_hist_{key}.png")
            plt.close(fig)
    elif plot_context == "qc":
        print("WARNING: skipping QC plots because total_counts is not present in adata.obs.")

    configure_plot_style()
    if should_plot_preprocess:
        save_preprocess_plots(adata, figures_dir)
        configure_plot_style()
    elif plot_context == "preprocess":
        print(
            "WARNING: skipping preprocess plots because preprocess outputs are missing "
            "(expected highly_variable, pca, or X_pca)."
        )

    if not should_plot_umap:
        return

    clustering_key = config.get("clustering", {}).get("key_added", "leiden")
    if "X_umap" not in adata.obsm:
        logging.warning("Skipping UMAP plots because X_umap is not present in adata.obsm.")
        return

    debug_max_cells = config.get("debug", {}).get("max_cells")
    umap_size = 3 if debug_max_cells is not None else 1
    palette = getattr(sc.pl.palettes, "godsnot_102", None)
    qc_umap_keys = [
        key
        for key in ["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo"]
        if key in adata.obs
    ]
    if "pct_counts_hb" in adata.obs:
        qc_umap_keys.append("pct_counts_hb")
    if qc_umap_keys:
        sc.pl.umap(
            adata,
            color=qc_umap_keys,
            size=umap_size,
            frameon=False,
            show=False,
            save="_qc_metrics.png",
        )

    overview_color = [
        key
        for key in [clustering_key, "sample_id", "condition", "batch"]
        if key in adata.obs
    ]
    if overview_color:
        sc.pl.umap(
            adata,
            color=overview_color,
            size=umap_size,
            frameon=False,
            show=False,
            save="_clusters_samples_conditions.png",
        )

    if clustering_key in adata.obs:
        sc.pl.umap(
            adata,
            color=clustering_key,
            palette=palette,
            figsize=(10, 8),
            size=umap_size,
            legend_loc="right margin",
            legend_fontoutline=2,
            add_outline=True,
            frameon=False,
            show=False,
            save="_leiden.png",
        )
        sc.pl.umap(
            adata,
            color=clustering_key,
            palette=palette,
            figsize=(10, 8),
            size=umap_size,
            legend_loc="on data",
            legend_fontsize=10,
            legend_fontoutline=2,
            add_outline=True,
            frameon=False,
            show=False,
            save="_leiden_labels.png",
        )

    if "individualID" in adata.obs:
        individual_categories = adata.obs["individualID"].astype("category").cat.categories
        individual_kwargs: dict[str, Any] = {
            "size": max(umap_size / 2, 1),
            "frameon": False,
            "show": False,
            "save": "_individualID.png",
        }
        if len(individual_categories) <= 30:
            individual_kwargs["legend_loc"] = "right margin"
        else:
            individual_kwargs["legend_loc"] = None
        sc.pl.umap(adata, color="individualID", **individual_kwargs)

    if "major.celltype" in adata.obs:
        sc.pl.umap(
            adata,
            color="major.celltype",
            palette=palette,
            size=umap_size,
            legend_loc="right margin",
            frameon=False,
            show=False,
            save="_major_celltype.png",
        )


def set_analysis_params(
    adata: ad.AnnData,
    config: dict[str, Any],
    config_path: Path,
    sample_sheet: Path | None,
    run_id: str | None = None,
    step: str | None = None,
) -> None:
    adata.uns["analysis_params"] = {
        "config": config,
        "config_path": str(config_path),
        "input_format": input_format(config),
        "sample_sheet": str(sample_sheet) if sample_sheet is not None else None,
        "run_id": run_id,
        "step": step,
    }


def load_latest_plot_input(paths: dict[str, Path]) -> ad.AnnData:
    for step in ("markers", "cluster", "harmony", "preprocess", "qc"):
        path = intermediate_h5ad_path(paths, step)
        if path.exists():
            logging.info("Loading latest available plotting input from %s", path)
            return sc.read_h5ad(path)
    raise FileNotFoundError(
        f"No step intermediates found in {paths['processed_dir']}. "
        "Run at least the 'cluster' step before plotting UMAPs."
    )


def load_plot_input(paths: dict[str, Path], plot_context: str) -> ad.AnnData:
    if plot_context == "auto":
        return load_latest_plot_input(paths)

    step_by_context = {
        "qc": "qc",
        "preprocess": "preprocess",
        "harmony": "harmony",
        "cluster": "cluster",
        "markers": "markers",
    }
    step = step_by_context[plot_context]
    path = intermediate_h5ad_path(paths, step)
    if not path.exists():
        raise FileNotFoundError(
            f"Missing plot input for context '{plot_context}': {path}. "
            f"Run the '{step}' step first."
        )
    logging.info("Loading %s plotting input from %s", plot_context, path)
    return sc.read_h5ad(path)


def load_cluster_input(paths: dict[str, Path]) -> tuple[ad.AnnData, Path]:
    harmony_path = intermediate_h5ad_path(paths, "harmony")
    if harmony_path.exists():
        logging.info("Loading clustering input from Harmony output: %s", harmony_path)
        return sc.read_h5ad(harmony_path), harmony_path

    preprocess_path = intermediate_h5ad_path(paths, "preprocess")
    if preprocess_path.exists():
        logging.info("Loading clustering input from preprocess output: %s", preprocess_path)
        return sc.read_h5ad(preprocess_path), preprocess_path

    raise FileNotFoundError(
        f"Missing clustering input. Expected {harmony_path} or {preprocess_path}. "
        "Run the 'harmony' or 'preprocess' step first."
    )


def run_step(args: argparse.Namespace, config: dict[str, Any]) -> None:
    run_id = resolve_run_id(config, args.step)
    paths = ensure_run_output_dirs(config, run_id)
    logging.info("Running step '%s' in %s", args.step, paths["run_dir"])
    if args.step != "plots" and args.plot_context != "auto":
        logging.warning(
            "Ignoring --plot-context=%s because --step is %s, not plots.",
            args.plot_context,
            args.step,
        )

    if args.step == "qc":
        adata, sample_sheet = load_input(args.config, config, args.sample_sheet)
        adata = apply_debug_subset(adata, config)
        add_qc_annotations(adata, config)
        adata, qc_summary = apply_filtering(adata, config)
        set_analysis_params(adata, config, args.config, sample_sheet, run_id, args.step)
        qc_summary.to_csv(paths["tables_dir"] / "qc_filtering_summary.csv", index=False)
        write_intermediate(adata, paths, args.step)
        return

    if args.step == "preprocess":
        adata = read_intermediate(paths, "qc")
        adata = preprocess(adata, config)
        set_analysis_params(adata, config, args.config, None, run_id, args.step)
        write_intermediate(adata, paths, args.step)
        return

    if args.step == "harmony":
        adata = read_intermediate(paths, "preprocess")
        harmony_summary = run_harmony_if_enabled(adata, config)
        set_analysis_params(adata, config, args.config, None, run_id, args.step)
        harmony_summary.to_csv(paths["tables_dir"] / "harmony_summary.csv", index=False)
        write_intermediate(adata, paths, args.step)
        return

    if args.step == "cluster":
        adata, input_h5ad = load_cluster_input(paths)
        clustering_summary = compute_embedding_and_clusters(adata, config, input_h5ad)
        set_analysis_params(adata, config, args.config, None, run_id, args.step)
        clustering_summary.to_csv(paths["tables_dir"] / "clustering_summary.csv", index=False)
        write_intermediate(adata, paths, args.step)
        return

    if args.step == "markers":
        marker_input = intermediate_h5ad_path(paths, "cluster")
        adata = read_intermediate(paths, "cluster")
        marker_table_path = paths["tables_dir"] / "cluster_markers.csv"
        marker_table, markers_summary = compute_markers(
            adata,
            config,
            marker_input,
            marker_table_path,
        )
        set_analysis_params(adata, config, args.config, None, run_id, args.step)
        marker_table.to_csv(marker_table_path, index=False)
        markers_summary.to_csv(paths["tables_dir"] / "markers_summary.csv", index=False)
        write_intermediate(adata, paths, args.step)
        return

    if args.step == "plots":
        adata = load_plot_input(paths, args.plot_context)
        save_plots(adata, paths, config, args.plot_context)
        set_analysis_params(adata, config, args.config, None, run_id, args.step)
        write_intermediate(adata, paths, args.step)
        return

    raise ValueError(f"Unsupported step: {args.step}")


def output_h5ad_path(paths: dict[str, Path], config: dict[str, Any]) -> Path:
    return paths["processed_dir"] / "processed_with_harmony.h5ad"


def main() -> None:
    setup_logging()
    args = parse_args()

    config = read_yaml(args.config)
    if args.step is not None:
        run_step(args, config)
        return
    if args.plot_context != "auto":
        logging.warning(
            "Ignoring --plot-context=%s because --step plots was not requested.",
            args.plot_context,
        )

    paths = ensure_output_dirs(config)

    adata, sample_sheet = load_input(args.config, config, args.sample_sheet)
    adata = apply_debug_subset(adata, config)
    add_qc_annotations(adata, config)
    adata, qc_summary = apply_filtering(adata, config)
    adata = preprocess(adata, config)
    harmony_summary = run_harmony_if_enabled(adata, config)
    clustering_summary = compute_embedding_and_clusters(adata, config)
    marker_table_path = paths["tables_dir"] / "cluster_markers.csv"
    marker_table, markers_summary = compute_markers(
        adata,
        config,
        None,
        marker_table_path,
    )

    set_analysis_params(adata, config, args.config, sample_sheet)

    qc_summary.to_csv(paths["tables_dir"] / "qc_filtering_summary.csv", index=False)
    harmony_summary.to_csv(paths["tables_dir"] / "harmony_summary.csv", index=False)
    clustering_summary.to_csv(paths["tables_dir"] / "clustering_summary.csv", index=False)
    marker_table.to_csv(marker_table_path, index=False)
    markers_summary.to_csv(paths["tables_dir"] / "markers_summary.csv", index=False)
    save_plots(adata, paths, config)
    adata.write_h5ad(output_h5ad_path(paths, config))


if __name__ == "__main__":
    main()

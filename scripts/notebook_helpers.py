"""Notebook-only helper functions for the interactive scRNA-seq workflow."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import subprocess
from typing import Optional

from IPython.display import display
import pandas as pd
import yaml


@dataclass
class NotebookContext:
    dataset_id: str
    project_root: Path
    python: Path
    data_root: Path
    data_path: Path
    run_id: str
    config_path: Path
    script_path: Path
    results_root: Path
    use_harmony: bool
    harmony_batch_key: str
    hvg_n_top_genes: int
    hvg_flavor: str
    neighbors_n_neighbors: int
    neighbors_n_pcs: int
    leiden_resolution: float
    random_state: int
    debug_max_cells: Optional[int]


_context: Optional[NotebookContext] = None


def configure_context(context: NotebookContext) -> None:
    global _context
    _context = context


def _export_helpers_to_namespace(notebook_globals: dict) -> None:
    helper_names = [
        "run_step",
        "write_current_config",
        "run_step_with_current_config",
        "get_latest_run",
        "show_workflow_status",
        "show_full_run_readiness_checklist",
        "suggest_next_action",
        "show_stage_decision_options",
        "show_qc_summary",
        "show_qc_figures",
        "show_qc_plots",
        "show_preprocess_summary",
        "show_preprocess_figures",
        "show_preprocess_plots",
        "show_harmony_summary",
        "show_harmony_figures",
        "show_harmony_plots",
        "show_clustering_summary",
        "show_clustering_figures",
        "show_clustering_plots",
        "show_marker_summary",
        "show_marker_table",
        "show_marker_figures",
        "show_marker_plots",
        "get_default_canonical_marker_panels",
        "show_canonical_marker_panels",
        "create_canonical_marker_review_request",
        "create_review_packet",
        "create_review_request",
        "create_review_bundle",
        "print_reviewer_prompt",
    ]
    module_globals = globals()
    for name in helper_names:
        notebook_globals[name] = module_globals[name]


def initialize_notebook_context(notebook_globals: dict) -> None:
    defaults = {
        "use_harmony": True,
        "harmony_batch_key": "individualID",
        "hvg_n_top_genes": 3000,
        "hvg_flavor": "seurat",
        "neighbors_n_neighbors": 15,
        "neighbors_n_pcs": 30,
        "leiden_resolution": 0.5,
        "random_state": 42,
        "debug_max_cells": 5000,
    }
    for name, value in defaults.items():
        notebook_globals.setdefault(name, value)

    context = NotebookContext(
        dataset_id=notebook_globals["DATASET_ID"],
        project_root=Path(notebook_globals["PROJECT_ROOT"]),
        python=Path(notebook_globals["PYTHON"]),
        data_root=Path(notebook_globals["DATA_ROOT"]),
        data_path=Path(notebook_globals["DATA_PATH"]),
        run_id=notebook_globals["RUN_ID"],
        config_path=Path(notebook_globals["CONFIG_PATH"]),
        script_path=Path(notebook_globals["SCRIPT_PATH"]),
        results_root=Path(notebook_globals["RESULTS_ROOT"]),
        use_harmony=notebook_globals["use_harmony"],
        harmony_batch_key=notebook_globals["harmony_batch_key"],
        hvg_n_top_genes=notebook_globals["hvg_n_top_genes"],
        hvg_flavor=notebook_globals["hvg_flavor"],
        neighbors_n_neighbors=notebook_globals["neighbors_n_neighbors"],
        neighbors_n_pcs=notebook_globals["neighbors_n_pcs"],
        leiden_resolution=notebook_globals["leiden_resolution"],
        random_state=notebook_globals["random_state"],
        debug_max_cells=notebook_globals["debug_max_cells"],
    )
    configure_context(context)
    _export_helpers_to_namespace(notebook_globals)
    print(f"Notebook helpers initialized for RUN_ID={context.run_id}")


DEFAULT_CANONICAL_MARKER_PANELS = {
    "Oligodendrocyte": ["MBP", "PLP1", "MOG", "MOBP", "MAG"],
    "OPC": ["PDGFRA", "CSPG4", "VCAN", "PTPRZ1"],
    "Astrocyte": ["AQP4", "GFAP", "ALDH1L1", "SLC1A3", "SLC1A2"],
    "Excitatory neuron": ["SLC17A7", "SLC17A6", "CAMK2A", "SATB2"],
    "Inhibitory neuron": ["GAD1", "GAD2", "SLC6A1", "DLX1", "DLX2"],
    "Microglia": ["CX3CR1", "P2RY12", "TMEM119", "AIF1", "CSF1R"],
    "Macrophage / immune": ["PTPRC", "LST1", "TYROBP", "C1QA", "C1QB", "MRC1"],
    "Endothelial": ["CLDN5", "FLT1", "PECAM1", "VWF"],
    "Pericyte / vascular mural": ["PDGFRB", "RGS5", "CSPG4", "ACTA2"],
    "Fibroblast / VLMC / ECM-like": ["COL1A1", "COL1A2", "DCN", "LUM", "NID1", "LAMA2"],
    "Ependymal / epithelial-like": ["FOXJ1", "PIFO", "TTR", "KRT18"],
    "Stress / immediate early": ["FOS", "JUN", "JUNB", "DUSP1"],
    "Hemoglobin / RBC warning": ["HBA1", "HBA2", "HBB"],
}


def get_default_canonical_marker_panels() -> dict[str, list[str]]:
    """Return editable canonical marker panels for broad brain cell-type review."""
    return {
        cell_class: list(markers)
        for cell_class, markers in DEFAULT_CANONICAL_MARKER_PANELS.items()
    }


def _canonical_marker_panel_frame(
    panels: Optional[dict[str, list[str]]] = None,
) -> pd.DataFrame:
    marker_panels = panels or get_default_canonical_marker_panels()
    return pd.DataFrame(
        [
            {
                "cell_class": cell_class,
                "canonical_markers": ", ".join(markers),
            }
            for cell_class, markers in marker_panels.items()
        ]
    )


def show_canonical_marker_panels(
    panels: Optional[dict[str, list[str]]] = None,
) -> None:
    display(_canonical_marker_panel_frame(panels))
    return None


def _ctx() -> NotebookContext:
    if _context is None:
        raise RuntimeError("Notebook helper context is not configured.")
    return _context


def check_paths(
    PROJECT_ROOT,
    PYTHON,
    DATA_ROOT,
    DATA_PATH,
    CONFIG_PATH,
    SCRIPT_PATH,
    RESULTS_ROOT,
):
    paths = {
        "PROJECT_ROOT": Path(PROJECT_ROOT),
        "PYTHON": Path(PYTHON),
        "DATA_ROOT": Path(DATA_ROOT),
        "DATA_PATH": Path(DATA_PATH),
        "CONFIG_PATH": Path(CONFIG_PATH),
        "SCRIPT_PATH": Path(SCRIPT_PATH),
        "RESULTS_ROOT": Path(RESULTS_ROOT),
    }
    return pd.DataFrame(
        [
            {
                "variable": variable,
                "path": str(path),
                "exists": path.exists(),
            }
            for variable, path in paths.items()
        ]
    )


def write_config_from_notebook(notebook_globals: dict):
    DATASET_ID = notebook_globals["DATASET_ID"]
    DATA_PATH = notebook_globals["DATA_PATH"]
    RUN_ID = notebook_globals["RUN_ID"]
    RESULTS_ROOT = notebook_globals["RESULTS_ROOT"]
    CONFIG_PATH = notebook_globals["CONFIG_PATH"]

    random_state = notebook_globals.get("random_state", 42)
    debug_max_cells = notebook_globals.get("debug_max_cells", None)
    normalization_target_sum = notebook_globals.get("normalization_target_sum", 10000)
    normalization_log1p = notebook_globals.get("normalization_log1p", True)
    scaling_max_value = notebook_globals.get("scaling_max_value", 10)
    hvg_n_top_genes = notebook_globals.get("hvg_n_top_genes", 3000)
    hvg_flavor = notebook_globals.get("hvg_flavor", "seurat")
    hvg_batch_key = notebook_globals.get("hvg_batch_key", None)
    pca_n_comps = notebook_globals.get("pca_n_comps", 50)
    pca_svd_solver = notebook_globals.get("pca_svd_solver", "arpack")
    use_harmony = notebook_globals.get("use_harmony", True)
    harmony_batch_key = notebook_globals.get("harmony_batch_key", "individualID")
    harmony_basis = notebook_globals.get("harmony_basis", "X_pca")
    harmony_adjusted_basis = notebook_globals.get(
        "harmony_adjusted_basis", "X_pca_harmony"
    )
    neighbors_n_neighbors = notebook_globals.get("neighbors_n_neighbors", 15)
    neighbors_n_pcs = notebook_globals.get("neighbors_n_pcs", 30)
    umap_random_state = notebook_globals.get("umap_random_state", random_state)
    leiden_resolution = notebook_globals.get("leiden_resolution", 0.5)
    clustering_key_added = notebook_globals.get("clustering_key_added", "leiden")
    clustering_metadata_keys = notebook_globals.get(
        "clustering_metadata_keys",
        ["leiden", "individualID", "major.celltype", "hcelltype", "region"],
    )
    markers_groupby = notebook_globals.get("markers_groupby", clustering_key_added)
    markers_method = notebook_globals.get("markers_method", "wilcoxon")
    markers_n_genes = notebook_globals.get("markers_n_genes", 100)
    marker_plot_top_n = notebook_globals.get("marker_plot_top_n", 5)
    qc_min_genes = notebook_globals.get("qc_min_genes", 200)
    qc_min_cells = notebook_globals.get("qc_min_cells", 10)
    qc_max_pct_mt = notebook_globals.get("qc_max_pct_mt", 8)
    qc_max_genes = notebook_globals.get("qc_max_genes", None)
    qc_max_counts = notebook_globals.get("qc_max_counts", None)
    qc_mt_gene_prefix = notebook_globals.get("qc_mt_gene_prefix", "MT-")
    qc_ribo_gene_prefixes = notebook_globals.get(
        "qc_ribo_gene_prefixes", ["RPS", "RPL"]
    )
    qc_hb_gene_pattern = notebook_globals.get("qc_hb_gene_pattern", "^HB[^(P)]")
    qc_mad_filter_enabled = notebook_globals.get("qc_mad_filter_enabled", True)
    qc_log1p_total_counts_nmads = notebook_globals.get(
        "qc_log1p_total_counts_nmads", 5
    )
    qc_log1p_n_genes_by_counts_nmads = notebook_globals.get(
        "qc_log1p_n_genes_by_counts_nmads", 5
    )
    qc_pct_counts_mt_nmads = notebook_globals.get("qc_pct_counts_mt_nmads", 3)
    qc_mt_direction = notebook_globals.get("qc_mt_direction", "upper")

    config = {
        "project_name": f"scrna_{DATASET_ID.lower()}",
        "input": {
            "format": "h5ad",
            "path": str(DATA_PATH),
        },
        "output": {
            "processed_dir": "data/processed",
            "figures_dir": "reports/figures",
            "tables_dir": "reports/tables",
            "logs_dir": "logs",
        },
        "run": {
            "id": RUN_ID,
            "add_timestamp": True,
            "base_dir": str(RESULTS_ROOT.resolve()),
        },
        "qc": {
            "min_genes": qc_min_genes,
            "min_cells": qc_min_cells,
            "max_pct_mt": qc_max_pct_mt,
            "max_genes": qc_max_genes,
            "max_counts": qc_max_counts,
            "mt_gene_prefix": qc_mt_gene_prefix,
            "ribo_gene_prefixes": qc_ribo_gene_prefixes,
            "hb_gene_pattern": qc_hb_gene_pattern,
            "mad_filter": {
                "enabled": qc_mad_filter_enabled,
                "log1p_total_counts_nmads": qc_log1p_total_counts_nmads,
                "log1p_n_genes_by_counts_nmads": qc_log1p_n_genes_by_counts_nmads,
                "pct_counts_mt_nmads": qc_pct_counts_mt_nmads,
                "mt_direction": qc_mt_direction,
            },
        },
        "normalization": {
            "target_sum": normalization_target_sum,
            "log1p": normalization_log1p,
        },
        "scaling": {
            "max_value": scaling_max_value,
        },
        "hvg": {
            "n_top_genes": hvg_n_top_genes,
            "flavor": hvg_flavor,
            "batch_key": hvg_batch_key,
        },
        "pca": {
            "n_comps": pca_n_comps,
            "svd_solver": pca_svd_solver,
            "random_state": random_state,
        },
        "harmony": {
            "enabled": use_harmony,
            "batch_key": harmony_batch_key,
            "basis": harmony_basis,
            "adjusted_basis": harmony_adjusted_basis,
        },
        "neighbors": {
            "n_neighbors": neighbors_n_neighbors,
            "n_pcs": neighbors_n_pcs,
        },
        "umap": {
            "random_state": umap_random_state,
        },
        "clustering": {
            "method": "leiden",
            "resolution": leiden_resolution,
            "key_added": clustering_key_added,
            "metadata_keys": clustering_metadata_keys,
            "random_state": random_state,
        },
        "debug": {
            "max_cells": debug_max_cells,
            "random_state": random_state,
        },
        "markers": {
            "groupby": markers_groupby,
            "method": markers_method,
            "n_genes": markers_n_genes,
            "plot_top_n": marker_plot_top_n,
        },
    }

    CONFIG_PATH.parent.mkdir(parents=True, exist_ok=True)
    with CONFIG_PATH.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(config, handle, sort_keys=False)
    return CONFIG_PATH


def write_current_config(notebook_globals: dict):
    """Write the server config from current notebook variables and refresh context."""
    config_path = write_config_from_notebook(notebook_globals)
    initialize_notebook_context(notebook_globals)
    print(f"Wrote current parameters to {config_path}")
    return config_path


def run_step(step, context=None):
    ctx = _ctx()
    command = [
        str(ctx.python),
        str(ctx.script_path),
        "--config",
        str(ctx.config_path),
        "--step",
        step,
    ]
    if context is not None:
        command.extend(["--plot-context", context])
    subprocess.run(command, check=True)


def run_step_with_current_config(notebook_globals: dict, step, context=None):
    """Write current notebook parameters, then run one pipeline step."""
    write_current_config(notebook_globals)
    run_step(step, context=context)


def get_latest_run():
    ctx = _ctx()
    run_dirs = sorted(
        [path for path in ctx.results_root.glob(f"{ctx.run_id}_*") if path.is_dir()],
        key=lambda path: path.name,
    )
    if not run_dirs:
        raise FileNotFoundError(
            f"No run folders found for {ctx.run_id} under {ctx.results_root}"
        )
    return run_dirs[-1]


def show_workflow_status():
    try:
        latest_run = get_latest_run()
    except FileNotFoundError as exc:
        print(f"No run folder found yet: {exc}")
        return None

    state_path = latest_run / "workflow_state.json"
    if state_path.exists():
        with state_path.open("r", encoding="utf-8") as handle:
            state = json.load(handle)
        status = pd.DataFrame(
            [
                {
                    "current_step": state.get("current_step"),
                    "completed_steps": ", ".join(state.get("completed_steps", [])),
                    "latest_h5ad": state.get("latest_h5ad"),
                    "available_tables": len(state.get("available_tables", [])),
                    "available_figures": len(state.get("available_figures", [])),
                }
            ]
        )
        print(f"Latest run: {latest_run}")
        display(status)
        return None

    expected_paths = [
        latest_run / "data" / "qc.h5ad",
        latest_run / "data" / "preprocessed.h5ad",
        latest_run / "data" / "harmony.h5ad",
        latest_run / "data" / "clustered.h5ad",
        latest_run / "data" / "markers.h5ad",
        latest_run / "tables" / "qc_filtering_summary.csv",
        latest_run / "tables" / "harmony_summary.csv",
        latest_run / "tables" / "clustering_summary.csv",
        latest_run / "tables" / "markers_summary.csv",
    ]
    status = pd.DataFrame(
        [
            {
                "artifact": path.name,
                "path": str(path.relative_to(latest_run)),
                "exists": path.exists(),
            }
            for path in expected_paths
        ]
    )
    print(f"Latest run: {latest_run}")
    display(status)
    return None


def _review_bundle_exists(latest_run: Optional[Path], step: str) -> tuple[str, str]:
    if latest_run is None:
        return "WARNING", "No run folder found yet."
    ctx = _ctx()
    bundle_dir = ctx.project_root / "prompts" / "review_requests" / f"{latest_run.name}_{step}"
    prompt_path = bundle_dir / "review_prompt.md"
    packet_path = bundle_dir / "biology_review_packet.md"
    if prompt_path.exists() and packet_path.exists():
        return "OK", f"{latest_run.name}_{step}"
    return "WARNING", f"Missing review bundle: prompts/review_requests/{latest_run.name}_{step}"


def _has_dangerous_harmony_key(batch_key: object) -> bool:
    key = str(batch_key or "").lower()
    dangerous_tokens = [
        "diagnosis",
        "disease",
        "condition",
        "phenotype",
        "outcome",
        "treatment",
        "cognition",
        "cognitive",
        "braak",
        "cerad",
        "pathologic",
        "pathological",
        "clinical",
        "case",
        "control",
    ]
    return any(token in key for token in dangerous_tokens)


def show_full_run_readiness_checklist():
    """Display non-executing checks before the user manually launches a full server run."""
    config = _read_current_config()
    try:
        latest_run = get_latest_run()
    except FileNotFoundError:
        latest_run = None

    checks = []

    debug_max_cells = (config.get("debug") or {}).get("max_cells")
    checks.append(
        {
            "check": "Debug mode disabled",
            "status": "OK" if debug_max_cells is None else "WARNING",
            "detail": (
                "debug.max_cells is None."
                if debug_max_cells is None
                else f"debug.max_cells is {debug_max_cells}; set it to None before the full run."
            ),
        }
    )

    checks.append(
        {
            "check": "Server config exists",
            "status": "OK" if _ctx().config_path.exists() else "WARNING",
            "detail": "configs/pipeline.server.yaml" if _ctx().config_path.exists() else "Missing configs/pipeline.server.yaml",
        }
    )

    for label, step in [
        ("QC reviewed", "qc"),
        ("Preprocess reviewed", "preprocess"),
        ("Harmony reviewed", "harmony"),
        ("Clustering reviewed", "cluster"),
        ("Markers reviewed", "markers"),
        ("Canonical marker review done", "canonical_markers"),
    ]:
        status, detail = _review_bundle_exists(latest_run, step)
        checks.append({"check": label, "status": status, "detail": detail})

    harmony = config.get("harmony") or {}
    harmony_enabled = bool(harmony.get("enabled", False))
    harmony_batch_key = harmony.get("batch_key")
    if harmony_enabled and _has_dangerous_harmony_key(harmony_batch_key):
        harmony_status = "WARNING"
        harmony_detail = (
            f"harmony.batch_key is {harmony_batch_key}; do not use biological, clinical, "
            "phenotype, diagnosis, condition, or outcome variables."
        )
    elif harmony_enabled:
        harmony_status = "OK"
        harmony_detail = f"Harmony enabled with batch key: {harmony_batch_key}"
    else:
        harmony_status = "INFO"
        harmony_detail = "Harmony disabled."
    checks.append(
        {
            "check": "Harmony batch key safety",
            "status": harmony_status,
            "detail": harmony_detail,
        }
    )

    latest_run_detail = latest_run.name if latest_run is not None else "No run folder found."
    checks.append({"check": "Latest reviewed run", "status": "INFO", "detail": latest_run_detail})

    display(pd.DataFrame(checks))

    key_parameters = {
        "run.id": (config.get("run") or {}).get("id"),
        "debug.max_cells": debug_max_cells,
        "qc.min_genes": (config.get("qc") or {}).get("min_genes"),
        "qc.max_pct_mt": (config.get("qc") or {}).get("max_pct_mt"),
        "hvg.n_top_genes": (config.get("hvg") or {}).get("n_top_genes"),
        "pca.n_comps": (config.get("pca") or {}).get("n_comps"),
        "harmony.enabled": harmony_enabled,
        "harmony.batch_key": harmony_batch_key,
        "neighbors.n_neighbors": (config.get("neighbors") or {}).get("n_neighbors"),
        "neighbors.n_pcs": (config.get("neighbors") or {}).get("n_pcs"),
        "clustering.resolution": (config.get("clustering") or {}).get("resolution"),
        "markers.groupby": (config.get("markers") or {}).get("groupby"),
    }
    display(
        pd.DataFrame(
            [{"parameter": key, "value": value} for key, value in key_parameters.items()]
        )
    )
    print("This checklist does not run analysis. Launch the full server run manually only after warnings are resolved.")
    return None


def suggest_next_action():
    try:
        latest_run = get_latest_run()
    except FileNotFoundError as exc:
        print(f"No run folder found yet: {exc}")
        return None

    state_path = latest_run / "workflow_state.json"
    if state_path.exists():
        with state_path.open("r", encoding="utf-8") as handle:
            state = json.load(handle)
        completed_stages = state.get("completed_steps", [])
        latest_completed = completed_stages[-1] if completed_stages else "none"
        completed = set(completed_stages)
    else:
        stage_paths = {
            "qc": latest_run / "data" / "qc.h5ad",
            "preprocess": latest_run / "data" / "preprocessed.h5ad",
            "harmony": latest_run / "data" / "harmony.h5ad",
            "cluster": latest_run / "data" / "clustered.h5ad",
            "markers": latest_run / "data" / "markers.h5ad",
        }
        completed_stages = [
            stage for stage, path in stage_paths.items() if path.exists()
        ]
        latest_completed = completed_stages[-1] if completed_stages else "none"
        completed = set(completed_stages)

    if "markers" in completed:
        next_stage = "review workflow / compare runs / annotation review"
        suggested_command = (
            "Use review workflow cells; compare runs or start annotation review."
        )
    elif "cluster" in completed:
        next_stage = "markers"
        suggested_command = 'run_step("markers")'
    elif "harmony" in completed:
        next_stage = "cluster"
        suggested_command = 'run_step("cluster")'
    elif "preprocess" in completed:
        next_stage = "harmony or clustering"
        suggested_command = 'run_step("harmony") or run_step("cluster")'
    elif "qc" in completed:
        next_stage = "preprocess"
        suggested_command = 'run_step("preprocess")'
    else:
        next_stage = "qc"
        suggested_command = 'run_step("qc")'

    print(f"Latest run: {latest_run}")
    print(f"Latest completed stage: {latest_completed}")
    print(f"Available next stage: {next_stage}")
    print(f"Suggested command: {suggested_command}")
    return None


def show_stage_decision_options(step):
    step_key = str(step).lower()
    if step_key == "qc":
        print("After QC review, choose one:")
        print("- Continue to preprocess if QC metrics and plots look acceptable.")
        print("- Edit QC parameters and rerun QC if filtering looks too strict or too permissive.")
        print("- Inspect additional QC diagnostics if high-count, high-gene, doublet, or per-sample balance concerns remain.")
        print("- Pause for human expert review if the QC decision is uncertain.")
        print("")
        print("Use suggest_next_action() only after making the QC decision.")
        return None
    if step_key == "preprocess":
        print("After preprocess review, choose one:")
        print("- Continue to Harmony if batch correction is justified and a valid technical batch key is available.")
        print("- Skip Harmony and proceed to clustering if integration is not justified for this run.")
        print("- Edit preprocess parameters and rerun preprocess if HVG selection, scaling, or PCA structure looks poor.")
        print("- Inspect additional PCA or metadata-colored diagnostics if technical effects remain unclear.")
        print("")
        print("Use suggest_next_action() only after making the preprocess decision.")
        return None
    if step_key == "harmony":
        print("After Harmony review, choose one:")
        print("- Continue to clustering with Harmony if the batch key is technical and integration looks justified.")
        print("- Skip Harmony and cluster on X_pca if integration is unnecessary or overcorrection risk is high.")
        print("- Edit Harmony parameters and rerun Harmony if the batch key or basis was wrong.")
        print("- Inspect additional batch/metadata plots if undercorrection or overcorrection remains unclear.")
        print("")
        print("Use suggest_next_action() only after making the Harmony decision.")
        return None
    if step_key in {"cluster", "clustering"}:
        print("After clustering review, choose one:")
        print("- Continue to marker analysis if UMAP structure and Leiden granularity look acceptable.")
        print("- Edit Leiden resolution and rerun clustering if clusters are too coarse or too granular.")
        print("- Edit neighbors_n_pcs or neighbors_n_neighbors and rerun clustering if graph structure looks unstable.")
        print("- Revisit Harmony if clusters appear batch-driven or if major biology appears overcorrected.")
        print("- Inspect additional metadata or QC overlays if cluster drivers remain unclear.")
        print("")
        print("Use suggest_next_action() only after making the clustering decision.")
        return None
    if step_key == "markers":
        print("After marker review, choose one:")
        print("- Continue to manual annotation only if marker genes are coherent and cluster quality looks acceptable.")
        print("- Recheck small, mixed, or low-quality clusters before assigning biological labels.")
        print("- Lower Leiden resolution and rerun clustering/markers if many clusters lack distinct markers.")
        print("- Increase Leiden resolution and rerun clustering/markers if broad clusters contain mixed marker programs.")
        print("- Revisit QC, Harmony, or clustering if marker patterns suggest doublets, stress, batch effects, or overcorrection.")
        print("")
        print("Use suggest_next_action() only after making the marker decision.")
        return None
    if step_key in {"canonical_markers", "canonical_marker_review", "annotation"}:
        print("After canonical marker review, choose one:")
        print("- Draft broad manual labels only for clusters with coherent marker evidence.")
        print("- Keep ambiguous, mixed, tiny, or low-confidence clusters unlabeled until expert review.")
        print("- Edit the canonical marker panel and rerun the review if important lineage markers are missing.")
        print("- Revisit clustering or Harmony if canonical markers suggest fragmentation, mixing, or overcorrection.")
        print("- Do not make disease, phenotype, or condition interpretation from this review.")
        return None

    print(f"No stage-specific decision options are configured for: {step}")
    return None


def show_qc_summary():
    latest_run = get_latest_run()
    qc_summary_path = latest_run / "tables" / "qc_filtering_summary.csv"
    if not qc_summary_path.exists():
        raise FileNotFoundError(f"QC summary not found: {qc_summary_path}")
    qc_summary = pd.read_csv(qc_summary_path)
    display(qc_summary)
    return None


def show_qc_figures():
    latest_run = get_latest_run()
    figures_dir = latest_run / "figures"
    expected_paths = [
        figures_dir / "violin_qc.png",
        figures_dir / "scatter_qc_counts_vs_mt.png",
        figures_dir / "qc_scatter_counts_vs_genes_mt.png",
    ]
    expected_paths.extend(
        sorted(figures_dir.glob("qc_hist*.png")) if figures_dir.exists() else []
    )

    print(f"Latest run: {latest_run}")
    print("QC figure paths:")
    for path in expected_paths:
        status = "FOUND" if path.exists() else "MISSING"
        print(f"- [{status}] {path}")


def show_qc_plots():
    from IPython.display import Image, display

    latest_run = get_latest_run()
    figures_dir = latest_run / "figures"
    qc_plots = [
        ("QC violin plot", figures_dir / "violin_qc.png"),
        ("Counts vs mitochondrial percentage", figures_dir / "scatter_qc_counts_vs_mt.png"),
        (
            "Counts vs genes colored by mitochondrial percentage",
            figures_dir / "qc_scatter_counts_vs_genes_mt.png",
        ),
        ("Mitochondrial percentage histogram", figures_dir / "qc_hist_pct_counts_mt.png"),
        ("Genes detected histogram", figures_dir / "qc_hist_n_genes_by_counts.png"),
        ("Total counts histogram", figures_dir / "qc_hist_total_counts.png"),
    ]

    for title, path in qc_plots:
        print(f"\n=== {title} ===")
        if path.exists():
            display(Image(filename=str(path)))
        else:
            print(f"WARNING: missing QC figure: {path}")


def show_preprocess_summary():
    import scanpy as sc

    latest_run = get_latest_run()
    preprocessed_path = latest_run / "data" / "preprocessed.h5ad"
    if not preprocessed_path.exists():
        raise FileNotFoundError(f"Preprocessed AnnData not found: {preprocessed_path}")
    adata = sc.read_h5ad(preprocessed_path)
    hvg_count = (
        int(adata.var["highly_variable"].sum())
        if "highly_variable" in adata.var
        else None
    )
    x_pca = adata.obsm.get("X_pca")
    summary = pd.DataFrame(
        [
            {
                "cells": adata.n_obs,
                "genes": adata.n_vars,
                "highly_variable_genes": hvg_count,
                "obsm_keys": ", ".join(adata.obsm.keys()),
                "has_X_pca": x_pca is not None,
                "X_pca_shape": str(x_pca.shape) if x_pca is not None else "",
            }
        ]
    )
    display(summary)
    return None


def show_preprocess_figures():
    latest_run = get_latest_run()
    figures_dir = latest_run / "figures"
    expected_paths = [
        figures_dir / "hvg_plot.png",
        figures_dir / "pca_scree.png",
        figures_dir / "pca_scatter_qc_metrics.png",
    ]

    print(f"Latest run: {latest_run}")
    print("Preprocess figure paths:")
    for path in expected_paths:
        status = "FOUND" if path.exists() else "MISSING"
        print(f"- [{status}] {path}")


def show_preprocess_plots():
    from IPython.display import Image, display

    latest_run = get_latest_run()
    figures_dir = latest_run / "figures"
    preprocess_plots = [
        ("Highly variable genes", figures_dir / "hvg_plot.png"),
        ("PCA variance explained", figures_dir / "pca_scree.png"),
        ("PCA colored by QC metrics", figures_dir / "pca_scatter_qc_metrics.png"),
    ]

    for title, path in preprocess_plots:
        print(f"\n=== {title} ===")
        if path.exists():
            display(Image(filename=str(path)))
        else:
            print(f"WARNING: missing preprocess figure: {path}")


def show_harmony_summary():
    latest_run = get_latest_run()
    summary_path = latest_run / "tables" / "harmony_summary.csv"
    if not summary_path.exists():
        raise FileNotFoundError(f"Harmony summary not found: {summary_path}")
    summary = pd.read_csv(summary_path)
    display(summary)
    return None


def _harmony_plot_paths():
    latest_run = get_latest_run()
    figures_dir = latest_run / "figures"
    if not figures_dir.exists():
        return []
    tokens = ["harmony", "integration", "batch"]
    return sorted(
        path
        for path in figures_dir.glob("*.png")
        if any(token in path.name.lower() for token in tokens)
    )


def show_harmony_figures():
    latest_run = get_latest_run()
    plot_paths = _harmony_plot_paths()
    print(f"Latest run: {latest_run}")
    if not plot_paths:
        print("No Harmony-related plot files found. Harmony plotting may be added later.")
        return None
    print("Harmony figure paths:")
    for path in plot_paths:
        print(f"- {path}")
    return None


def show_harmony_plots():
    from IPython.display import Image, display

    plot_paths = _harmony_plot_paths()
    if not plot_paths:
        print("No Harmony plots are currently generated.")
        return None
    for path in plot_paths:
        print(f"\n=== {path.name} ===")
        display(Image(filename=str(path)))
    return None


def show_clustering_summary():
    latest_run = get_latest_run()
    summary_path = latest_run / "tables" / "clustering_summary.csv"
    if not summary_path.exists():
        raise FileNotFoundError(f"Clustering summary not found: {summary_path}")
    summary = pd.read_csv(summary_path)
    display(summary)
    return None


def show_clustering_figures():
    latest_run = get_latest_run()
    figures_dir = latest_run / "figures"
    expected_paths = [
        figures_dir / "umap_qc_counts_genes.png",
        figures_dir / "umap_leiden.png",
        figures_dir / "umap_leiden_numbered.png",
        figures_dir / "umap_metadata_panel.png",
    ]
    print(f"Latest run: {latest_run}")
    print("Clustering figure paths:")
    for path in expected_paths:
        status = "FOUND" if path.exists() else "MISSING"
        print(f"- [{status}] {path}")
    return None


def show_clustering_plots():
    from IPython.display import Image, display

    latest_run = get_latest_run()
    figures_dir = latest_run / "figures"
    clustering_plots = [
        ("UMAP colored by counts and genes", figures_dir / "umap_qc_counts_genes.png"),
        ("UMAP colored by Leiden clusters", figures_dir / "umap_leiden.png"),
        ("UMAP with Leiden labels on data", figures_dir / "umap_leiden_numbered.png"),
        ("UMAP metadata sanity-check panel", figures_dir / "umap_metadata_panel.png"),
    ]
    for title, path in clustering_plots:
        print(f"\n=== {title} ===")
        if path.exists():
            display(Image(filename=str(path)))
        else:
            print(f"WARNING: optional clustering figure is missing: {path}")
    return None


def show_marker_summary():
    latest_run = get_latest_run()
    summary_path = latest_run / "tables" / "markers_summary.csv"
    if not summary_path.exists():
        raise FileNotFoundError(f"Marker summary not found: {summary_path}")
    summary = pd.read_csv(summary_path)
    display(summary)
    return None


def show_marker_table(n=50):
    latest_run = get_latest_run()
    marker_table_path = latest_run / "tables" / "cluster_markers.csv"
    if not marker_table_path.exists():
        raise FileNotFoundError(f"Marker table not found: {marker_table_path}")
    marker_table = pd.read_csv(marker_table_path)
    display(marker_table.head(n))
    return None


def show_marker_figures():
    latest_run = get_latest_run()
    figures_dir = latest_run / "figures"
    expected_paths = [
        figures_dir / "marker_dotplot_top_genes.png",
        figures_dir / "marker_matrixplot_top_genes.png",
    ]
    print(f"Latest run: {latest_run}")
    print("Marker figure paths:")
    for path in expected_paths:
        status = "FOUND" if path.exists() else "MISSING"
        print(f"- [{status}] {path}")
    return None


def show_marker_plots():
    from IPython.display import Image, display

    latest_run = get_latest_run()
    figures_dir = latest_run / "figures"
    marker_plots = [
        ("Marker dotplot top genes", figures_dir / "marker_dotplot_top_genes.png"),
        ("Marker matrixplot top genes", figures_dir / "marker_matrixplot_top_genes.png"),
    ]
    for title, path in marker_plots:
        print(f"\n=== {title} ===")
        if path.exists():
            display(Image(filename=str(path)))
        else:
            print(f"WARNING: optional marker figure is missing: {path}")
    return None


def _read_current_config() -> dict:
    ctx = _ctx()
    if not ctx.config_path.exists():
        return {}
    with ctx.config_path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def _append_parameter_section(lines: list[str], title: str, values: dict) -> None:
    if not values:
        return
    lines.append(f"### {title}")
    for key, value in values.items():
        lines.append(f"- {key}: `{value}`")
    lines.append("")


def _flatten_parameters(section: dict, prefix: str = "") -> dict:
    flattened = {}
    for key, value in (section or {}).items():
        label = f"{prefix}.{key}" if prefix else str(key)
        if isinstance(value, dict):
            flattened.update(_flatten_parameters(value, label))
        else:
                flattened[label] = value
    return flattened


def _safe_review_table(path: Path, drop_columns: Optional[set[str]] = None) -> pd.DataFrame:
    table = pd.read_csv(path)
    columns_to_drop = set(drop_columns or set())
    columns_to_drop.update(
        column
        for column in table.columns
        if column.lower().endswith("_path") or column.lower() in {"input_h5ad"}
    )
    return table.drop(columns=[column for column in columns_to_drop if column in table.columns])


def _stage_parameter_lines(step: str) -> list[str]:
    config = _read_current_config()
    step_key = str(step).lower()
    lines = ["## Current Stage Parameters", ""]

    if step_key == "qc":
        _append_parameter_section(lines, "QC", _flatten_parameters(config.get("qc", {})))
        _append_parameter_section(lines, "Debug", _flatten_parameters(config.get("debug", {})))
    elif step_key == "preprocess":
        _append_parameter_section(
            lines, "Normalization", _flatten_parameters(config.get("normalization", {}))
        )
        _append_parameter_section(lines, "Scaling", _flatten_parameters(config.get("scaling", {})))
        _append_parameter_section(lines, "Highly Variable Genes", _flatten_parameters(config.get("hvg", {})))
        _append_parameter_section(lines, "PCA", _flatten_parameters(config.get("pca", {})))
        _append_parameter_section(lines, "Debug", _flatten_parameters(config.get("debug", {})))
    elif step_key == "harmony":
        _append_parameter_section(lines, "Harmony", _flatten_parameters(config.get("harmony", {})))
    elif step_key in {"cluster", "clustering"}:
        _append_parameter_section(lines, "Neighbors", _flatten_parameters(config.get("neighbors", {})))
        _append_parameter_section(lines, "UMAP", _flatten_parameters(config.get("umap", {})))
        _append_parameter_section(lines, "Clustering", _flatten_parameters(config.get("clustering", {})))
    elif step_key in {"markers", "canonical_markers", "canonical_marker_review", "annotation"}:
        _append_parameter_section(lines, "Markers", _flatten_parameters(config.get("markers", {})))
        _append_parameter_section(lines, "Clustering Context", _flatten_parameters(config.get("clustering", {})))
    else:
        lines.append("Stage-specific parameters not configured for this review step.")
        lines.append("")

    if lines == ["## Current Stage Parameters", ""]:
        lines.append("No current config parameters were found for this review step.")
        lines.append("")

    return lines


def _stage_figure_paths(step: str, figures_dir: Path) -> list[Path]:
    all_figure_paths = sorted(figures_dir.glob("*.png")) if figures_dir.exists() else []
    step_key = str(step).lower()
    if step_key == "qc":
        return [
            path
            for path in all_figure_paths
            if path.name
            in {
                "violin_qc.png",
                "scatter_qc_counts_vs_mt.png",
                "qc_scatter_counts_vs_genes_mt.png",
            }
            or path.name.startswith("qc_hist_")
        ]
    if step_key == "preprocess":
        return [
            path
            for path in all_figure_paths
            if path.name in {"hvg_plot.png", "pca_scree.png", "pca_scatter_qc_metrics.png"}
        ]
    if step_key == "harmony":
        return [
            path
            for path in all_figure_paths
            if any(token in path.name.lower() for token in ["harmony", "integration", "batch"])
        ]
    if step_key in {"cluster", "clustering"}:
        return [
            path
            for path in all_figure_paths
            if path.name
            in {
                "umap_qc_counts_genes.png",
                "umap_leiden.png",
                "umap_leiden_numbered.png",
                "umap_metadata_panel.png",
            }
        ]
    if step_key == "markers":
        return [
            path
            for path in all_figure_paths
            if path.name
            in {
                "marker_dotplot_top_genes.png",
                "marker_matrixplot_top_genes.png",
            }
        ]
    if step_key in {"canonical_markers", "canonical_marker_review", "annotation"}:
        return [
            path
            for path in all_figure_paths
            if path.name
            in {
                "marker_dotplot_top_genes.png",
                "marker_matrixplot_top_genes.png",
            }
        ]
    return all_figure_paths


def _best_practice_checklist_text(step: str) -> str:
    ctx = _ctx()
    step_key = str(step).lower()
    checklist_by_step = {
        "qc": "qc.md",
        "preprocess": "preprocess.md",
        "harmony": "harmony.md",
        "cluster": "clustering.md",
        "clustering": "clustering.md",
        "markers": "markers.md",
        "canonical_markers": "canonical_annotation.md",
        "canonical_marker_review": "canonical_annotation.md",
        "annotation": "canonical_annotation.md",
    }
    checklist_name = checklist_by_step.get(step_key)
    if checklist_name is None:
        return "Stage-specific best-practice checklist is not configured for this step."
    checklist_path = ctx.project_root / "prompts" / "best_practices" / checklist_name
    if not checklist_path.exists():
        return f"Stage-specific best-practice checklist not found: `{checklist_path.name}`"
    return checklist_path.read_text(encoding="utf-8").strip()


def _append_qc_doublet_summary(lines: list[str], obs: pd.DataFrame) -> None:
    doublet_columns = [
        column
        for column in ("is.doublet", "is_doublet", "doublet")
        if column in obs.columns
    ]
    score_columns = [
        column
        for column in ("doublet_score", "scrublet_score", "doublet.score")
        if column in obs.columns
    ]
    if not doublet_columns and not score_columns:
        return

    lines.append("## Doublet Summary")
    if doublet_columns:
        doublet_column = doublet_columns[0]
        doublet_values = obs[doublet_column]
        doublet_bool = doublet_values.astype(str).str.lower().isin(
            {"true", "1", "yes", "doublet"}
        )
        called_doublets = int(doublet_bool.sum())
        total_with_calls = int(doublet_values.notna().sum())
        doublet_pct = (
            round(called_doublets / total_with_calls * 100, 2)
            if total_with_calls
            else 0.0
        )
        lines.append(f"- doublet call column: `{doublet_column}`")
        lines.append(f"- called doublets: `{called_doublets}` / `{total_with_calls}` (`{doublet_pct}%`)")

        required_qc = {"total_counts", "n_genes_by_counts"}
        if required_qc.issubset(obs.columns) and total_with_calls:
            high_counts_cutoff = obs["total_counts"].quantile(0.95)
            high_genes_cutoff = obs["n_genes_by_counts"].quantile(0.95)
            high_tail = obs[
                (obs["total_counts"] >= high_counts_cutoff)
                | (obs["n_genes_by_counts"] >= high_genes_cutoff)
            ]
            if not high_tail.empty:
                high_tail_doublet_bool = doublet_bool.loc[high_tail.index]
                high_tail_doublets = int(high_tail_doublet_bool.sum())
                high_tail_pct = round(high_tail_doublets / len(high_tail) * 100, 2)
                lines.append(
                    "- high-count/high-gene tail doublets "
                    f"(top 5% by either metric): `{high_tail_doublets}` / "
                    f"`{len(high_tail)}` (`{high_tail_pct}%`)"
                )
    if score_columns:
        score_column = score_columns[0]
        score = pd.to_numeric(obs[score_column], errors="coerce").dropna()
        if not score.empty:
            lines.append(f"- doublet score column: `{score_column}`")
            lines.append(
                "- doublet score summary: "
                f"median `{score.median():.4g}`, "
                f"95th percentile `{score.quantile(0.95):.4g}`, "
                f"max `{score.max():.4g}`"
            )
    lines.append("")


def create_review_packet(
    step,
    contact_sheet_name: Optional[str] = None,
    canonical_marker_panels: Optional[dict[str, list[str]]] = None,
):
    ctx = _ctx()
    latest_run = get_latest_run()
    packet_path = latest_run / "biology_review_packet.md"
    tables_dir = latest_run / "tables"
    figures_dir = latest_run / "figures"
    data_dir = latest_run / "data"

    qc_summary_path = tables_dir / "qc_filtering_summary.csv"
    harmony_summary_path = tables_dir / "harmony_summary.csv"
    clustering_summary_path = tables_dir / "clustering_summary.csv"
    markers_summary_path = tables_dir / "markers_summary.csv"
    marker_table_path = tables_dir / "cluster_markers.csv"
    h5ad_paths = sorted(data_dir.glob("*.h5ad")) if data_dir.exists() else []
    latest_h5ad = h5ad_paths[-1] if h5ad_paths else None
    step_key = str(step).lower()
    figure_paths = _stage_figure_paths(step_key, figures_dir)

    expected_h5ad_by_step = {
        "qc": data_dir / "qc.h5ad",
        "preprocess": data_dir / "preprocessed.h5ad",
        "harmony": data_dir / "harmony.h5ad",
        "cluster": data_dir / "clustered.h5ad",
        "clustering": data_dir / "clustered.h5ad",
        "markers": data_dir / "markers.h5ad",
        "canonical_markers": data_dir / "markers.h5ad",
        "canonical_marker_review": data_dir / "markers.h5ad",
        "annotation": data_dir / "markers.h5ad",
        "plots": data_dir / "plots.h5ad",
    }
    expected_h5ad = expected_h5ad_by_step.get(step_key)
    selected_h5ad = (
        expected_h5ad
        if expected_h5ad is not None and expected_h5ad.exists()
        else latest_h5ad
    )
    marker_review_steps = {"markers", "canonical_markers", "canonical_marker_review", "annotation"}
    if step_key in marker_review_steps and expected_h5ad is not None and not expected_h5ad.exists():
        clustered_h5ad = data_dir / "clustered.h5ad"
        selected_h5ad = clustered_h5ad if clustered_h5ad.exists() else latest_h5ad
    h5ad_warning = None
    if expected_h5ad is not None and not expected_h5ad.exists():
        if step_key in marker_review_steps:
            h5ad_warning = (
                f"Expected h5ad for step '{step}' was missing: `{expected_h5ad}`. "
                "Falling back to clustered h5ad if available."
            )
        else:
            h5ad_warning = (
                f"Expected h5ad for step '{step}' was missing: `{expected_h5ad}`. "
                "Falling back to latest available h5ad."
            )

    lines = []
    lines.append("# Biology Review Packet")
    lines.append("")
    lines.append(f"- Review step: `{step}`")
    lines.append(f"- Run: `{latest_run.name}`")
    lines.append("- Config: `configs/pipeline.server.yaml`")
    lines.append("")
    lines.extend(_stage_parameter_lines(step))

    lines.append("## QC Summary")
    if qc_summary_path.exists():
        lines.append(pd.read_csv(qc_summary_path).to_markdown(index=False))
    else:
        lines.append("QC summary not available.")
    lines.append("")

    if step_key == "harmony":
        lines.append("## Harmony Summary")
        if harmony_summary_path.exists():
            lines.append(_safe_review_table(harmony_summary_path).to_markdown(index=False))
        else:
            lines.append("Harmony summary not available.")
        lines.append("")

    if step_key in {"cluster", "clustering"}:
        lines.append("## Clustering Summary")
        if clustering_summary_path.exists():
            lines.append(_safe_review_table(clustering_summary_path).to_markdown(index=False))
        else:
            lines.append("Clustering summary not available.")
        lines.append("")

    if step_key in marker_review_steps:
        lines.append("## Marker Summary")
        if markers_summary_path.exists():
            lines.append(_safe_review_table(markers_summary_path).to_markdown(index=False))
        else:
            lines.append("Marker summary not available.")
        lines.append("")

    if step_key in {"canonical_markers", "canonical_marker_review", "annotation"}:
        lines.append("## Canonical Marker Panels")
        lines.append(_canonical_marker_panel_frame(canonical_marker_panels).to_markdown(index=False))
        lines.append("")
        lines.append("Interpretation rules:")
        lines.append("- Use canonical panels as supporting evidence, not automatic labels.")
        lines.append("- Prefer broad labels unless marker evidence is strong and specific.")
        lines.append("- Flag mixed, tiny, low-quality, or conflicting clusters as uncertain.")
        lines.append("- Do not infer disease, condition, phenotype, or treatment effects.")
        lines.append("")

    lines.append("## Latest AnnData Summary")
    if h5ad_warning is not None:
        lines.append(f"WARNING: {h5ad_warning}")
    if selected_h5ad is not None:
        try:
            import scanpy as sc

            adata = sc.read_h5ad(selected_h5ad)
            lines.append(f"- Latest h5ad: `{selected_h5ad.name}`")
            lines.append(f"- Cells: `{adata.n_obs}`")
            lines.append(f"- Genes: `{adata.n_vars}`")
            obs_columns = list(adata.obs.columns)
            shown_obs_columns = obs_columns[:25]
            lines.append(
                f"- obs columns (showing first {len(shown_obs_columns)} of {len(obs_columns)}): "
                + ", ".join(f"`{col}`" for col in shown_obs_columns)
            )
            if step_key in {"preprocess", "harmony", "cluster", "clustering", "markers", "canonical_markers", "canonical_marker_review", "annotation"}:
                hvg_count = (
                    int(adata.var["highly_variable"].sum())
                    if "highly_variable" in adata.var
                    else "not available"
                )
                x_pca = adata.obsm.get("X_pca")
                lines.append(f"- highly variable genes: `{hvg_count}`")
                lines.append(f"- raw present: `{adata.raw is not None}`")
                lines.append(f"- X_pca present: `{x_pca is not None}`")
                lines.append(
                    f"- X_pca shape: `{x_pca.shape if x_pca is not None else 'not available'}`"
                )
                if step_key in {"harmony", "cluster", "clustering", "markers", "canonical_markers", "canonical_marker_review", "annotation"}:
                    x_pca_harmony = adata.obsm.get("X_pca_harmony")
                    lines.append(f"- X_pca_harmony present: `{x_pca_harmony is not None}`")
                    lines.append(
                        f"- X_pca_harmony shape: `{x_pca_harmony.shape if x_pca_harmony is not None else 'not available'}`"
                    )
                if step_key in {"cluster", "clustering", "markers", "canonical_markers", "canonical_marker_review", "annotation"}:
                    x_umap = adata.obsm.get("X_umap")
                    clustering_key = "leiden"
                    lines.append(f"- X_umap present: `{x_umap is not None}`")
                    lines.append(
                        f"- X_umap shape: `{x_umap.shape if x_umap is not None else 'not available'}`"
                    )
                    lines.append(f"- Leiden key present: `{clustering_key in adata.obs}`")
                    if clustering_key in adata.obs:
                        cluster_sizes = adata.obs[clustering_key].astype(str).value_counts().sort_index()
                        preview = "; ".join(
                            f"{cluster}:{count}" for cluster, count in cluster_sizes.head(10).items()
                        )
                        if len(cluster_sizes) > 10:
                            preview += f"; ... ({len(cluster_sizes)} total)"
                        lines.append(f"- number of clusters: `{len(cluster_sizes)}`")
                        lines.append(f"- cluster sizes preview: `{preview}`")
            lines.append("- obsm keys: " + ", ".join(f"`{key}`" for key in adata.obsm.keys()))
            lines.append("- layers: " + ", ".join(f"`{key}`" for key in adata.layers.keys()))
            if step_key == "qc":
                _append_qc_doublet_summary(lines, adata.obs)
        except Exception as exc:
            lines.append(f"Latest h5ad found, but summary failed: `{exc}`")
    else:
        lines.append("Latest h5ad not available.")
    lines.append("")

    if step_key in marker_review_steps:
        lines.append("## Marker Table")
        if marker_table_path.exists():
            marker_table = pd.read_csv(marker_table_path)
            if "cluster" in marker_table.columns:
                preview = marker_table.groupby("cluster", sort=False).head(5)
            else:
                preview = marker_table.head(50)
            lines.append(preview.to_markdown(index=False))
        else:
            lines.append("Marker table not available for this step.")
        lines.append("")

    lines.append("## Figure Summary")
    if contact_sheet_name is not None:
        lines.append(f"Review the attached contact sheet: `{contact_sheet_name}`")
    elif figure_paths:
        lines.append("Stage figure files are available in the pipeline results folder.")
    else:
        lines.append("Figure files not available for this step.")
    lines.append("")
    lines.append("This packet intentionally excludes raw expression matrix data.")

    packet_path.write_text("\n".join(lines), encoding="utf-8")
    return packet_path


def _build_reviewer_prompt(
    step, packet_text: str, contact_sheet_name: Optional[str] = None
) -> str:
    step_guidance = ""
    if step == "qc":
        step_guidance = """
## QC Review Guidance

Please evaluate:

- Whether QC filtering appears overly aggressive or too permissive.
- Whether mitochondrial filtering looks biologically reasonable.
- Whether QC distributions suggest technical artifacts.
- Whether the retained cell count looks reasonable.
- Whether preprocessing can safely proceed.
""".strip()
    elif step == "preprocess":
        step_guidance = """
## Preprocess Review Guidance

Please evaluate:

- Whether normalization and log1p look appropriate.
- Whether HVG selection looks reasonable.
- Whether PCA variance structure looks defensible.
- Whether technical effects remain visible in PCA.
- Whether Harmony or clustering can safely proceed.
""".strip()
    elif step == "harmony":
        step_guidance = """
## Harmony Review Guidance

Please evaluate:

- Whether Harmony was justified for this dataset.
- Whether the selected batch key is appropriate.
- Whether the batch key might encode biological or clinical variation.
- Whether Harmony ran or was correctly skipped.
- Whether overcorrection risk should be checked before clustering.
- Whether clustering should use X_pca_harmony or original X_pca.
""".strip()
    elif step in {"cluster", "clustering"}:
        step_guidance = """
## Clustering Review Guidance

Please evaluate:

- Whether UMAP structure looks plausible.
- Whether Leiden resolution appears too coarse or too granular.
- Whether clusters appear driven by QC metrics.
- Whether Harmony may have overcorrected or undercorrected.
- Whether neighbor graph parameters should be adjusted.
- Whether Leiden resolution should be increased or decreased.
- Whether marker analysis can safely proceed.
""".strip()
    elif step == "markers":
        step_guidance = """
## Marker Review Guidance

Please evaluate:

- Whether marker genes look biologically plausible.
- Whether clusters appear well-separated by markers.
- Whether any clusters look low-quality, doublet-like, or mixed.
- Whether Leiden resolution should be adjusted.
- Whether manual annotation can be considered later.
- Whether any cluster should be merged, split, or rechecked.
""".strip()
    elif step in {"canonical_markers", "canonical_marker_review", "annotation"}:
        step_guidance = """
## Canonical Marker Review Guidance

Please evaluate:

- Which broad canonical cell class is supported for each cluster, if any.
- Which canonical markers support each tentative label.
- Which clusters are ambiguous, mixed, tiny, low-quality, or potentially doublet-like.
- Whether any clusters should remain unlabeled pending human expert review.
- Whether clustering, Harmony, or marker parameters should be revisited before annotation.

Do not make disease, phenotype, condition, treatment, diagnosis, or outcome interpretations.
Do not present labels as final; use tentative broad labels with confidence notes.
""".strip()

    contact_sheet_guidance = ""
    if contact_sheet_name is not None:
        contact_sheet_guidance = (
            "\nPlease inspect the attached contact sheet image: "
            f"`{contact_sheet_name}`.\n"
        )

    best_practice_checklist = _best_practice_checklist_text(step)

    prompt = f"""
Use prompts/biology_reviewer.md as your reviewer instructions.

Please review this single-cell run using only the packet content below and any plots/screenshots I provide. Do not access data/ or inspect h5ad files directly unless I explicitly allow it.
{contact_sheet_guidance}

{step_guidance}

## Stage-Specific Best-Practice Checklist

{best_practice_checklist}

{packet_text}
""".strip()

    return prompt


def _create_contact_sheet(figure_paths: list[Path], output_path: Path) -> Optional[Path]:
    if not figure_paths:
        return None

    import math
    import matplotlib.image as mpimg
    import matplotlib.pyplot as plt

    columns = min(2, len(figure_paths))
    rows = math.ceil(len(figure_paths) / columns)
    fig, axes = plt.subplots(rows, columns, figsize=(7 * columns, 5 * rows))
    try:
        axes = list(axes.ravel())
    except AttributeError:
        axes = [axes]

    for axis, path in zip(axes, figure_paths):
        image = mpimg.imread(path)
        axis.imshow(image)
        axis.set_title(path.name, fontsize=10)
        axis.axis("off")

    for axis in axes[len(figure_paths) :]:
        axis.axis("off")

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return output_path


def create_review_bundle(step, canonical_marker_panels: Optional[dict[str, list[str]]] = None):
    ctx = _ctx()
    latest_run = get_latest_run()
    step_key = str(step).lower()
    bundle_dir = ctx.project_root / "prompts" / "review_requests" / f"{latest_run.name}_{step_key}"
    bundle_dir.mkdir(parents=True, exist_ok=True)

    figure_paths = _stage_figure_paths(step_key, latest_run / "figures")
    contact_sheet_name = f"plots_contact_sheet_{latest_run.name}_{step_key}.png"
    contact_sheet_path = _create_contact_sheet(figure_paths, bundle_dir / contact_sheet_name)

    packet_path = create_review_packet(
        step,
        contact_sheet_name=contact_sheet_name if contact_sheet_path is not None else None,
        canonical_marker_panels=canonical_marker_panels,
    )
    packet_text = packet_path.read_text(encoding="utf-8")
    prompt = _build_reviewer_prompt(
        step,
        packet_text,
        contact_sheet_name=contact_sheet_name if contact_sheet_path is not None else None,
    )
    request_path = bundle_dir / "review_prompt.md"
    bundle_packet_path = bundle_dir / "biology_review_packet.md"
    request_path.write_text(prompt + "\n", encoding="utf-8")
    bundle_packet_path.write_text(packet_text + "\n", encoding="utf-8")

    try:
        relative_bundle_dir = bundle_dir.relative_to(ctx.project_root)
        relative_request_path = request_path.relative_to(ctx.project_root)
        relative_contact_sheet_path = (
            contact_sheet_path.relative_to(ctx.project_root)
            if contact_sheet_path is not None
            else None
        )
    except ValueError:
        relative_bundle_dir = bundle_dir
        relative_request_path = request_path
        relative_contact_sheet_path = contact_sheet_path

    print(f"Review bundle written to: {bundle_dir}")
    print(f"Stage plots included in contact sheet: {len(figure_paths)}")
    print("")
    print("After syncing prompts/review_requests locally, in Codex CLI send:")
    if relative_contact_sheet_path is not None:
        print(f"Review @{relative_request_path} @{relative_contact_sheet_path}")
    else:
        print(f"Review @{relative_request_path}")
    print("")
    print(f"Bundle folder: {relative_bundle_dir}")
    if contact_sheet_path is not None:
        from IPython.display import Image, display

        display(Image(filename=str(contact_sheet_path)))
    return bundle_dir


def create_review_request(step):
    return create_review_bundle(step)


def create_canonical_marker_review_request(
    panels: Optional[dict[str, list[str]]] = None,
):
    return create_review_bundle("canonical_markers", canonical_marker_panels=panels)


def print_reviewer_prompt(step):
    packet_path = create_review_packet(step)
    packet_text = packet_path.read_text(encoding="utf-8")
    prompt = _build_reviewer_prompt(step, packet_text)
    print(prompt)

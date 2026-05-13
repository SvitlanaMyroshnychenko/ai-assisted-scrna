# Notebook Usage Guide

This guide explains how to use `notebooks/01_run_scrna_pipeline.ipynb` as a controlled cockpit for the scRNA-seq workflow.

The notebook is not meant to be a black-box auto-analysis button. It is a human-in-the-loop workflow: run one stage, inspect outputs, ask Codex for review, decide whether to continue or change parameters, then move to the next stage.

## Big Picture

The project has three main parts:

- `scripts/scrna_pipeline.py`: deterministic pipeline engine.
- `scripts/notebook_helpers.py`: notebook orchestration helpers.
- `notebooks/01_run_scrna_pipeline.ipynb`: interactive workflow cockpit.

The intended workflow is:

1. Set run-wide parameters.
2. Write the initial config.
3. Run QC.
4. Review QC.
5. Decide whether to continue or change QC parameters.
6. Run preprocess.
7. Review preprocess.
8. Run Harmony if justified.
9. Review Harmony.
10. Run clustering.
11. Review clustering.
12. Run markers.
13. Review markers.
14. Run canonical marker review.
15. Use the full-run readiness checklist before a full dataset run.

## Important Safety Rules

- Do not store raw count matrices or real private metadata in this local project.
- Do not put credentials, private server paths, or patient identifiers into committed files.
- Do not run full analysis locally.
- The full analysis is run manually on the remote server by the user.
- Codex should not SSH into the server unless explicitly asked.
- Review prompts should contain summaries, plots, and marker tables, not raw expression matrices.

## Before Starting

Open the notebook:

```text
notebooks/01_run_scrna_pipeline.ipynb
```

Run the initialization cells at the top of the notebook.

If you changed `scripts/notebook_helpers.py` or deployed updated code to the server, restart the notebook kernel before continuing. Jupyter keeps old Python modules in memory, so parameter fixes may not take effect until restart.

## Section 1: Initialize Helpers

This section imports the notebook helper module and exposes helper functions as `nh.*`.

Run this section first.

If this section was not run, later cells such as `nh.run_step_with_current_config(...)` will fail.

## Section 2: Initial Parameters

This section defines run-wide parameters such as:

- dataset ID
- run ID
- data path
- results root
- Python executable
- config path
- random state
- debug cell limit

For debug runs, `debug_max_cells` is usually a number such as `5000`.

For a full run, `debug_max_cells` should be `None`.

Do not switch to full run until the full-run readiness checklist looks acceptable.

## Section 3: Write Config

This writes `configs/pipeline.server.yaml` from the current notebook parameters.

This is only a config snapshot. Stage-specific run cells also rewrite the config before launching their own stage.

The important idea:

- You can change parameters inside a stage section.
- The stage run cell writes those current values into the config.
- Then the pipeline stage runs using the updated config.

## Section 4: QC

QC is the first real analysis stage.

Typical parameters include:

- `qc_min_genes`
- `qc_min_cells`
- `qc_max_pct_mt`
- `qc_max_genes`
- `qc_max_counts`
- MAD filter settings
- mitochondrial/ribosomal/hemoglobin marker settings

Run order:

1. Edit QC parameters if needed.
2. Run the QC stage cell.
3. Run QC plots.
4. Inspect QC summary and plots.
5. Create QC review request.
6. Send the generated Codex CLI command.
7. Decide whether to continue or change QC parameters.

Do not continue just because the code ran. Continue only after QC looks technically reasonable.

Common reasons to rerun QC:

- too many cells removed
- too few cells removed
- high mitochondrial tail remains
- high-count/high-gene cells look suspicious
- doublet summary looks concerning
- per-sample balance looks questionable

## Section 5: Preprocess

Preprocess includes:

- normalization
- log transform
- highly variable gene selection
- scaling
- PCA

Typical parameters include:

- `normalization_target_sum`
- `normalization_log1p`
- `hvg_n_top_genes`
- `hvg_flavor`
- `hvg_batch_key`
- `scaling_max_value`
- `pca_n_comps`
- `pca_svd_solver`

Run order:

1. Edit preprocess parameters if needed.
2. Run preprocess.
3. Generate preprocess plots.
4. Review HVG and PCA diagnostics.
5. Create preprocess review request.
6. Decide whether Harmony or clustering should follow.

Important:

The number of PCA components computed in preprocess is not always the same as the number used in neighbors. For example, you may compute 55 PCs but use only 30 for neighbors.

## Section 6: Harmony

Harmony is optional and should only be used when batch correction is justified.

Typical parameters include:

- `use_harmony`
- `harmony_batch_key`
- `harmony_basis`
- `harmony_adjusted_basis`

Safe Harmony keys are technical batch variables.

Dangerous Harmony keys include variables related to:

- diagnosis
- disease
- condition
- phenotype
- outcome
- cognition
- treatment
- clinical state

Run order:

1. Decide whether Harmony is justified.
2. Set Harmony parameters.
3. Run Harmony.
4. Generate Harmony diagnostic plots.
5. Create Harmony review request.
6. Decide whether clustering should use Harmony or original PCA.

Harmony can help remove technical donor or batch effects, but it can also erase real biology. Always review it.

## Section 7: Clustering

Clustering includes:

- neighbor graph construction
- UMAP
- Leiden clustering

Typical parameters include:

- `neighbors_n_neighbors`
- `neighbors_n_pcs`
- `umap_random_state`
- `leiden_resolution`
- `clustering_key_added`
- `clustering_metadata_keys`

Run order:

1. Edit clustering parameters if needed.
2. Run clustering.
3. Generate clustering plots.
4. Review UMAP, cluster sizes, QC overlays, and metadata overlays.
5. Create clustering review request.
6. Decide whether marker analysis can proceed.

Common reasons to rerun clustering:

- clusters are too coarse
- clusters are too granular
- many tiny clusters appear
- UMAP is driven by QC metrics
- clusters appear batch-driven
- Harmony appears to overcorrect

## Section 8: Markers

Marker analysis tests marker genes for the current clusters.

Typical parameters include:

- `markers_groupby`
- `markers_method`
- `markers_n_genes`
- `marker_plot_top_n`

Run order:

1. Edit marker parameters if needed.
2. Run marker analysis.
3. Generate marker plots.
4. Inspect marker summary.
5. Inspect marker table.
6. Inspect dotplot and matrixplot.
7. Create marker review request.
8. Decide whether clusters are ready for canonical marker review.

Marker analysis does not automatically assign cell types.

The correct interpretation is:

- markers provide evidence
- Codex can suggest tentative interpretations
- the user decides what to trust
- final labels should not be written automatically without human review

## Section 9: Canonical Marker Review

Canonical marker review is an interpretation checkpoint after markers.

It does not run a new analysis step.

It does not write labels into AnnData.

It creates a review request asking Codex to compare cluster markers against broad canonical marker panels.

Default panels include broad brain cell classes such as:

- oligodendrocytes
- OPCs
- astrocytes
- excitatory neurons
- inhibitory neurons
- microglia
- macrophage/immune cells
- endothelial cells
- pericyte/vascular mural cells
- fibroblast/VLMC/ECM-like cells
- ependymal/epithelial-like cells
- stress markers
- hemoglobin/RBC warning markers

Run order:

1. Create or edit `canonical_marker_panels`.
2. Display the panels.
3. Create canonical marker review request.
4. Send the generated Codex CLI command.
5. Read Codex response.
6. Decide whether tentative broad labels are reasonable.

Important:

Canonical marker review should produce tentative broad labels and uncertainty notes. It should not produce final annotation automatically.

## Section 10: Full Run Readiness Checklist

This section checks whether the run is ready to move from debug mode to full run.

It does not launch the pipeline.

It checks:

- whether `debug.max_cells` is still enabled
- whether `configs/pipeline.server.yaml` exists
- whether review bundles exist for each stage
- whether Harmony batch key looks dangerous
- which latest run is being checked
- key final parameters

Before a full run:

- `debug.max_cells` should be `None`
- all major stage reviews should be completed
- Harmony decision should be explicit
- clustering parameters should be accepted
- marker review should not show major technical problems

## Review Requests and Codex CLI

After each stage, the notebook creates a review bundle in:

```text
prompts/review_requests/
```

Each bundle usually contains:

- `review_prompt.md`
- `biology_review_packet.md`
- `plots_contact_sheet_<run>_<stage>.png`

The notebook prints a Codex CLI command like:

```text
Review @prompts/review_requests/<run>_<stage>/review_prompt.md @prompts/review_requests/<run>_<stage>/plots_contact_sheet_<run>_<stage>.png
```

Run that command in Codex CLI after syncing the local `prompts/review_requests` folder.

The prompt is designed to avoid raw data and avoid full private server paths where possible.

## Contact Sheets

The pipeline writes individual plots into:

```text
results/<analysis_id>/figures/
```

The review bundle writes only a contact sheet into:

```text
prompts/review_requests/<analysis_id>_<stage>/
```

This keeps the review folder compact while still letting Codex CLI inspect plots.

## If Parameters Do Not Change

If you edit a parameter in the notebook but the run still uses the old value:

1. Make sure you edited the parameter in the correct stage section.
2. Make sure you ran the stage cell that uses `nh.run_step_with_current_config(...)`.
3. Restart the Jupyter kernel if helper code was changed or redeployed.
4. Rerun initialization cells.
5. Rerun the relevant stage.

This problem often happens because Jupyter kept an old version of `notebook_helpers.py` in memory.

## Debug Run vs Full Run

Use debug runs to test:

- whether parameters behave reasonably
- whether plots are generated
- whether review prompts are useful
- whether the workflow can proceed stage by stage

Use full runs only after:

- QC is accepted
- preprocess is accepted
- Harmony decision is accepted
- clustering is accepted
- marker review is acceptable
- full-run readiness checklist warnings are resolved

## What Not To Do

Do not skip reviews just because a stage completed.

Do not use disease, diagnosis, condition, phenotype, cognition, treatment, or outcome variables as Harmony batch keys.

Do not assign final cell types only from one marker table without reviewing canonical markers and cluster quality.

Do not interpret disease biology before technical validation and annotation are complete.

Do not run the full dataset while `debug.max_cells` is still set to a small number.

## Minimal Happy Path

For a simple debug iteration:

1. Run initialization.
2. Set initial parameters.
3. Write config.
4. Run QC.
5. Generate QC review request.
6. Accept or adjust QC.
7. Run preprocess.
8. Generate preprocess review request.
9. Decide Harmony.
10. Run Harmony if needed.
11. Generate Harmony review request.
12. Run clustering.
13. Generate clustering review request.
14. Run markers.
15. Generate marker review request.
16. Run canonical marker review.
17. Check full-run readiness.

The goal is not speed. The goal is controlled, reproducible analysis with human decisions at every risky point.

# Project Context

## Goal

Build a controlled, human-in-the-loop AI-assisted scRNA-seq/snRNA-seq workflow.

The workflow should help a human analyst run a reproducible Scanpy-based pipeline stage by stage, review each stage with Codex, adjust parameters when needed, and only then continue.

This is not an autonomous agent that runs analysis by itself. The notebook is a workflow cockpit, and the user remains responsible for decisions.

## Execution Model

- Codex runs locally on the user's computer.
- Sensitive biological data lives on the remote server.
- Full analysis is executed manually by the user on the remote server.
- The local project is code/templates/prompts/documentation only.
- Codex must not SSH into the server unless explicitly asked.
- Codex must not run analysis without explicit user permission.
- Local validation is limited to lightweight checks such as Python syntax and notebook JSON validity.

## Current Architecture

- `scripts/scrna_pipeline.py`: deterministic pipeline engine.
- `scripts/notebook_helpers.py`: notebook orchestration and review helpers.
- `notebooks/01_run_scrna_pipeline.ipynb`: interactive workflow cockpit.
- `configs/pipeline.template.yaml`: local template config.
- `configs/pipeline.server.yaml`: generated runtime config for server execution, should remain uncommitted if it contains real paths.
- `prompts/review_requests/`: local review bundles for Codex CLI.
- `NOTEBOOK_USAGE_GUIDE.md`: detailed user guide for the notebook workflow.

## Pipeline Stages

The deterministic pipeline supports these main stages:

1. QC
2. Preprocess
3. Harmony
4. Clustering
5. Markers
6. Plots

The notebook exposes these as a controlled staged workflow rather than a single monolithic run.

## Current Notebook Structure

1. Initialize Helpers
2. Initial Parameters
3. Write Config
4. QC
5. Preprocess
6. Harmony
7. Clustering
8. Markers
9. Canonical Marker Review
10. Full Run Readiness Checklist

The old empty `Review Workflow` section has been removed/replaced by the readiness checklist.

## Stage-Specific Parameter UX

Stage-specific parameters are intentionally placed near the stage where they matter.

The stage run cells use:

```python
nh.run_step_with_current_config(globals(), "<stage>")
```

This writes the current notebook values into `configs/pipeline.server.yaml` immediately before running that stage.

This pattern is used so the user can change parameters, rerun only the relevant stage, review the new result, and then decide whether to continue.

## Review Workflow

After each major stage, the notebook can create a local review bundle:

```text
prompts/review_requests/<run_folder>_<stage>/
```

Typical bundle contents:

- `review_prompt.md`
- `biology_review_packet.md`
- `plots_contact_sheet_<run_folder>_<stage>.png`

The notebook prints a Codex CLI command that references the prompt and contact sheet.

Review bundles are designed to include summaries, plots, and marker previews, not raw data.

Absolute server paths are avoided in generated review prompts where possible. Summary tables are sanitized before insertion into review packets so path-like columns such as `input_h5ad` and `*_path` do not leak full server paths.

## Contact Sheets

Pipeline plots remain in:

```text
results/<run_folder>/figures/
```

Review bundles contain only a contact sheet image, not duplicate copies of every individual plot.

This keeps review requests compact while still allowing Codex CLI to inspect stage plots.

## Current Implemented Review Points

### QC

QC has stage-local parameters and review outputs.

The review packet includes:

- QC filtering summary
- current QC parameters
- AnnData summary
- QC plots/contact sheet
- optional doublet summary when doublet columns exist

QC plots were polished with clearer histograms, scatter plots, and threshold overlays.

### Preprocess

Preprocess has stage-local parameters for:

- normalization
- log transform
- HVG selection
- scaling
- PCA

The review focuses on HVG behavior, PCA structure, QC metric influence, and whether Harmony/clustering can proceed.

### Harmony

Harmony has stage-local parameters for:

- enabled/disabled
- batch key
- input PCA basis
- adjusted PCA basis

Harmony diagnostic plots compare before/after PCA representations colored by batch key, QC metrics, and selected metadata where available.

The workflow explicitly warns against using biological, clinical, disease, phenotype, condition, cognition, treatment, or outcome variables as Harmony batch keys.

### Clustering

Clustering has stage-local parameters for:

- neighbors
- number of PCs for graph construction
- UMAP random state
- Leiden resolution
- clustering key
- metadata overlays

Clustering plots include UMAP cluster views, numbered clusters, QC overlays, and configured metadata panels.

### Markers

Marker analysis has stage-local parameters for:

- marker groupby key
- marker method
- number of ranked genes
- number of top genes shown in plots

Marker plots include:

- marker dotplot
- marker matrixplot

Marker review does not assign final cell types automatically. It evaluates whether marker evidence is coherent enough to proceed.

### Canonical Marker Review

Canonical marker review is an interpretation checkpoint after marker analysis.

It does not run a new analysis step.

It does not write labels into AnnData.

It compares cluster marker genes against editable broad canonical marker panels for brain cell classes and warning programs.

Default panels include:

- oligodendrocyte
- OPC
- astrocyte
- excitatory neuron
- inhibitory neuron
- microglia
- macrophage/immune
- endothelial
- pericyte/vascular mural
- fibroblast/VLMC/ECM-like
- ependymal/epithelial-like
- stress/immediate early
- hemoglobin/RBC warning

The expected output is tentative broad labels with uncertainty notes, not final automated annotation.

### Full Run Readiness Checklist

The final notebook section provides a non-executing checklist before moving from debug subset to full run.

It checks:

- whether debug mode is still enabled
- whether `configs/pipeline.server.yaml` exists
- whether review bundles exist for major stages
- whether the Harmony batch key looks risky
- latest reviewed run folder
- key final parameters

It does not launch the pipeline.

## Important Design Decisions

- Keep pipeline behavior deterministic and script-driven.
- Keep notebook as orchestration cockpit, not primary analysis logic.
- Keep helpers in `scripts/notebook_helpers.py`.
- Prefer small, safe refactors.
- Do not rewrite into LangChain, LangGraph, or a web app.
- Do not make the workflow autonomous.
- Human decisions gate progression between risky stages.
- Do not assign final cell types automatically.
- Do not interpret disease biology before technical validation and annotation review.

## Current Refactor Status

Completed:

- notebook helper extraction
- stage-specific config writing
- QC parameter handling
- preprocess parameter handling
- Harmony parameter handling
- clustering parameter handling
- marker parameter handling
- review bundle generation
- local review request folder workflow
- contact sheet generation
- review packet path sanitization
- QC plot improvements
- Harmony diagnostic plots
- clustering plot improvements
- marker dotplot/matrixplot improvements
- canonical marker review section
- full-run readiness checklist
- detailed notebook usage guide

Likely next areas:

- optional run comparison utilities
- optional parameter change log
- optional final run summary/report
- optional manual annotation draft table

## Safety Constraints

- Do not commit raw data.
- Do not commit `.h5ad`, `.h5`, `.loom`, `.mtx`, FASTQ, BAM, CRAM, or processed matrices.
- Do not commit real patient/sample identifiers or sensitive metadata.
- Do not commit credentials or private server paths.
- Do not inspect raw biological data unless explicitly allowed.
- Do not run full analysis locally.
- Do not SSH into the server unless explicitly asked.
- Do not run or modify server-side analysis without explicit user permission.
- Do not delete, overwrite, or move raw data.
- Do not use destructive commands unless the user explicitly approves the exact action.

## Definition of Done For Notebook Refactors

A notebook refactor is acceptable when:

- it preserves pipeline behavior;
- it keeps parameters editable at the appropriate stage;
- it does not require raw data locally;
- it does not trigger analysis during implementation;
- it does not add server access;
- generated review prompts avoid raw data and unnecessary absolute paths;
- lightweight local validation passes.

# Project: local Codex agent for remote scRNA-seq workflow

You are helping develop a reproducible, human-in-the-loop scRNA-seq/snRNA-seq workflow.

This project is a controlled analysis cockpit, not an autonomous biological analysis system. The human user decides when to run analysis, when to change parameters, and how to interpret biological results.

## Execution Model

- Codex runs locally on the user's computer.
- Sensitive scRNA-seq data lives only on the remote server.
- The local project is a code/template/prompt/documentation project and should not contain raw data.
- Full analysis is executed manually by the user on the remote server.
- Do not assume raw data exists locally.
- Do not attempt to SSH into the server unless the user explicitly asks.
- Do not run analysis unless the user explicitly asks.
- If the user says "do not change anything", "nothing do", "only explain", "do not edit", or similar, only analyze/explain and do not edit files.

## Required First Reads

When starting substantial work in this project, read:

- `AGENTS.md`
- `PROJECT_CONTEXT.md`
- `README.md`

For notebook usage details, also consult:

- `NOTEBOOK_USAGE_GUIDE.md`

## Current Architecture

- `scripts/scrna_pipeline.py`: deterministic Scanpy pipeline engine.
- `scripts/notebook_helpers.py`: notebook orchestration and review helpers.
- `notebooks/01_run_scrna_pipeline.ipynb`: interactive workflow cockpit.
- `configs/pipeline.template.yaml`: committed local template config.
- `configs/pipeline.server.yaml`: generated runtime config for server execution; should remain uncommitted if it contains real paths.
- `prompts/biology_reviewer.md`: general biology reviewer instructions.
- `prompts/review_requests/`: local review bundles for Codex CLI.
- `PROJECT_CONTEXT.md`: current architecture and refactor status.
- `NOTEBOOK_USAGE_GUIDE.md`: user-facing notebook guide.

## Preferred Stack

- Python
- Scanpy
- AnnData
- harmonypy
- pandas
- numpy
- scipy
- matplotlib
- PyYAML
- Jupyter notebooks for inspection
- Scripts for reproducible analysis

## Data Safety Rules

- Never request raw count matrices, `.h5ad`, `.h5`, `.loom`, `.mtx`, FASTQ, BAM, CRAM, or patient metadata.
- Never add real sample identifiers, patient identifiers, private paths, private server names, private project names, or access credentials to committed files.
- Do not open, print, summarize, or upload large biological data files.
- Do not inspect raw data, `data/`, generated `.h5ad`, or ignored result directories unless the user explicitly allows it.
- Respect `.codexignore`; never inspect ignored files or directories.
- Never run full analysis locally.
- Never write commands that delete, overwrite, move, or modify raw data.
- Never use destructive commands such as `rm -rf` unless the user explicitly approves the exact command.
- Do not commit `.env`, `.env.*`, logs, generated results, large matrices, processed `.h5ad` files, figures, or tables.

## Input and Config Rules

- The current notebook workflow uses `input.format: h5ad` with `input.path: /path/to/dataset.h5ad`.
- Sample sheets are optional and are only used for `input.format: 10x_mtx`.
- Use `configs/samples.template.csv` locally as a committed example for `10x_mtx` multi-sample input.
- The real server sample sheet for `10x_mtx` workflows must be named `configs/samples.server.csv` and must remain gitignored.
- `configs/samples.server.csv` is not required for the current `.h5ad` notebook workflow.
- Use `configs/pipeline.template.yaml` locally.
- The real/generated server pipeline config must be named `configs/pipeline.server.yaml` and must remain gitignored if it contains real server paths.

## Project Layout

Expected files and folders:

```text
AGENTS.md
README.md
PROJECT_CONTEXT.md
NOTEBOOK_USAGE_GUIDE.md
.gitignore
env/environment.yml
configs/pipeline.template.yaml
configs/samples.template.csv
scripts/scrna_pipeline.py
scripts/notebook_helpers.py
notebooks/01_run_scrna_pipeline.ipynb
prompts/biology_reviewer.md
prompts/review_requests/
reports/
results/
logs/
data/
```

## Notebook Workflow

The current notebook structure is:

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

The notebook is a workflow cockpit. Pipeline logic belongs in `scripts/scrna_pipeline.py`; orchestration helpers belong in `scripts/notebook_helpers.py`.

## Stage-Specific Parameter Rules

- Keep stage-specific parameters near the notebook section where they are used.
- Do not hard-code all downstream analysis parameters only at the beginning of the notebook.
- Stage run cells should write the current notebook parameters immediately before executing the stage.
- Prefer `nh.run_step_with_current_config(globals(), "<stage>")` for stage runs.
- A user should be able to edit parameters for a stage, rerun that stage, review the result, and decide whether to continue.
- Do not silently change pipeline behavior while refactoring notebook UX.

## Standard scRNA-seq Pipeline Order

1. Read parameters from YAML config.
2. Read sample metadata from CSV only for `input.format: 10x_mtx`.
3. Load data (`h5ad` or 10x mtx).
4. Add metadata to `adata.obs` when sample metadata is provided.
5. Concatenate samples for multi-sample 10x input.
6. Compute QC metrics.
7. Filter cells and genes.
8. Normalize.
9. Log transform.
10. Select highly variable genes.
11. Scale.
12. PCA.
13. Harmony if enabled.
14. Neighbors.
15. UMAP.
16. Leiden clustering.
17. Marker genes.
18. Save outputs.

## Best-Practices Reference

Use the following resource as a methodological reference:

- https://www.sc-best-practices.org/preamble.html
- Heumos et al., Nature Reviews Genetics (2023)

Guidelines:

- Treat this as guidance, not rigid law.
- Prefer best-practice workflows when possible.
- Do not hard-code thresholds blindly.
- Always justify QC thresholds.
- Avoid over-interpreting biological results.
- Warn about overcorrection in integration methods such as Harmony.
- Prefer stage-specific best-practice checklists over copying large external chapters into prompts.
- If uncertain, recommend additional diagnostics rather than confident conclusions.

## QC Rules

- Compute QC metrics before filtering.
- Record thresholds.
- Do not hard-code thresholds blindly.
- Make mitochondrial detection configurable.
- Save QC plots.
- Include threshold overlays and useful diagnostics where possible.
- Flag high-count/high-gene tails, mitochondrial concerns, doublet concerns, and per-sample imbalance as review items.

## Normalization and Preprocess Rules

- Preserve raw counts in `adata.layers["counts"]`.
- Use config-driven parameters.
- Avoid destructive overwrites.
- Keep normalization, HVG, scaling, and PCA parameters stage-local in the notebook.
- Review HVG behavior and PCA structure before clustering.

## Harmony Rules

- Run Harmony after PCA and before neighbors.
- Use only valid technical batch keys.
- Never use biological condition, disease, diagnosis, treatment, phenotype, cognition, clinical status, or outcome variables as Harmony batch keys.
- Save Harmony embedding in `adata.obsm["X_pca_harmony"]`.
- Use the adjusted embedding in neighbors only when Harmony is accepted.
- Warn about overcorrection.
- Encourage comparison of pre/post-Harmony structure and metadata overlays.

## Clustering and Marker Rules

- Parameters must come from config.
- Save clustering summaries and marker tables.
- Use review plots to evaluate cluster granularity, QC-driven structure, batch-driven structure, and small clusters.
- Marker analysis is evidence generation, not final annotation.
- Do not assign cell types without marker evidence and human review.

## Canonical Marker Review Rules

- Canonical marker review is an interpretation checkpoint after marker analysis.
- It should not run a new analysis step.
- It should not write labels into AnnData automatically.
- It should compare marker results with editable canonical marker panels.
- It should produce tentative broad labels, supporting evidence, uncertainty notes, and warnings for ambiguous/mixed/tiny clusters.
- Keep disease, condition, phenotype, treatment, diagnosis, and outcome interpretation out of canonical marker review.

## Review Bundle and Prompt Safety

- Review bundles live under `prompts/review_requests/`.
- Review packets may include summaries, marker table previews, and contact sheets.
- Review packets must not include raw expression matrices.
- Avoid full server paths in generated review prompts.
- Sanitize path-like columns such as `input_h5ad` and `*_path` before inserting summary tables into prompts.
- Prefer contact sheets over copying every individual plot into review request folders.
- The review folder should remain compact and safe to sync locally.

## Biology Reviewer Mode

- Use `prompts/biology_reviewer.md` when the user asks for a biology review of a run.
- Review only `biology_review_packet.md`, marker tables, QC summaries, and user-provided plots or screenshots unless the user explicitly allows deeper inspection.
- Do not access `data/` or inspect `.h5ad` files directly during biology review unless explicitly allowed.
- Provide tentative interpretations only, with marker evidence and uncertainty clearly separated.
- Recommend practical next checks such as continuing, changing Leiden resolution, changing the Harmony batch key, skipping Harmony, or inspecting specific clusters.
- Never recommend biological disease, diagnosis, treatment, phenotype, condition, or outcome variables as Harmony batch keys.
- Do not make disease or phenotype claims from clustering or marker review.

## Differential Expression

- Distinguish marker analysis from condition differential expression.
- Prefer replicate-aware or pseudobulk methods for condition DE.
- Do not treat cells as independent biological replicates.
- Do not recommend disease/condition DE from debug clustering alone.

## Full Run Readiness

- Full run should happen only after debug-stage review is acceptable.
- `debug.max_cells` should be `None` before a full dataset run.
- The readiness checklist is non-executing; it must not launch the pipeline.
- Do not add a one-click full-run button without explicit user request.
- The user manually launches full server analysis.

## Reproducibility

- Prefer scripts over notebook-only logic.
- Keep pipeline parameters config-driven.
- Save parameters in `adata.uns["analysis_params"]`.
- Save summary tables and figures.
- Keep review prompts and contact sheets reproducible from run outputs.

## Local Validation

Allowed:

```bash
python -m py_compile scripts/*.py
```

Also allowed:

- notebook JSON validity checks
- lightweight documentation checks

Do not run full analysis locally.

## Server Execution

Pipeline will run on the server manually:

```bash
conda run -n scrna-agent python scripts/scrna_pipeline.py --config configs/pipeline.server.yaml
```

Do not execute this command unless the user explicitly asks.

## Codex Behavior

- Make minimal, safe changes.
- Do not fabricate data.
- Do not run heavy computations.
- Ask before risky operations.
- If the user asks only for advice, do not edit files.
- If the user permits edits, stay within local project files unless explicitly told otherwise.
- Never touch the server without explicit permission.
- Preserve pipeline behavior unless the user asks for a behavioral change.
- Prefer small refactors that improve notebook UX, review safety, or reproducibility.

## Definition of Done

A change is acceptable when:

- pipeline behavior is preserved or the requested behavior change is explicit;
- stage parameters remain config-driven;
- notebook UX remains stage-local and review-driven;
- no raw data is accessed;
- no server access occurs;
- generated review prompts avoid raw data and unnecessary absolute paths;
- lightweight local validation passes when relevant;
- documentation is updated when workflow behavior changes.

# Project: local Codex agent for remote scRNA-seq workflow

You are helping develop a reproducible scRNA-seq analysis workflow.

## Execution model

- Codex runs locally on the user's computer.
- Sensitive scRNA-seq data lives only on the remote server.
- The local project is a code/template project and should not contain raw data.
- Full analysis is executed manually by the user on the remote server.
- Do not assume raw data exists locally.
- Do not attempt to SSH into the server unless the user explicitly asks.

## Preferred stack

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

## Data safety rules

- Never request raw count matrices, `.h5ad`, `.h5`, `.loom`, `.mtx`, FASTQ, BAM, CRAM, or patient metadata.
- Never add real sample identifiers, patient identifiers, private paths, private server names, private project names, or access credentials to committed files.
- The current notebook workflow uses `input.format: h5ad` with `input.path: /path/to/dataset.h5ad`.
- Sample sheets are optional and are only used for `input.format: 10x_mtx`.
- Use `configs/samples.template.csv` locally as a committed example for `10x_mtx` multi-sample input.
- The real server sample sheet for `10x_mtx` workflows must be named `configs/samples.server.csv` and must remain gitignored.
- `configs/samples.server.csv` is not required for the current `.h5ad` notebook workflow.
- Use `configs/pipeline.template.yaml` locally.
- The real server pipeline config must be named `configs/pipeline.server.yaml` and must remain gitignored.
- Never run full analysis locally.
- Never write commands that delete, overwrite, move, or modify raw data.
- Never use `rm -rf` unless the user explicitly approves the exact command.
- Do not open, print, summarize, or upload large biological data files.
- Do not commit `.env`, `.env.*`, logs, generated results, large matrices, or processed `.h5ad` files.
- Respect `.codexignore`; never inspect ignored files or directories
- 
## Project layout

Expected files and folders:

AGENTS.md  
README.md  
.gitignore  
env/environment.yml  
configs/pipeline.template.yaml  
configs/samples.template.csv  
scripts/scrna_pipeline.py 
notebooks/01_run_scrna_pipeline.ipynb  
reports/  
results/  
logs/  
data/

## Best-practices reference

Use the following resource as a methodological reference:

- https://www.sc-best-practices.org/preamble.html
- Heumos et al., Nature Reviews Genetics (2023)

Guidelines:

- Treat this as guidance, not rigid rules.
- Prefer best-practice workflows when possible.
- Do not hard-code thresholds blindly.
- Always justify QC thresholds.
- Avoid over-interpreting biological results.
- Warn about overcorrection in integration methods like Harmony.

## Standard scRNA-seq pipeline

1. Read parameters from YAML config
2. Read sample metadata from CSV only for `input.format: 10x_mtx`
3. Load data (`h5ad` or 10x mtx)
4. Add metadata to adata.obs when sample metadata is provided
5. Concatenate samples for multi-sample 10x input
6. Compute QC metrics
7. Filter cells and genes
8. Normalize
9. Log transform
10. Select highly variable genes
11. Scale
12. PCA
13. Harmony (if enabled)
14. Neighbors
15. UMAP
16. Leiden clustering
17. Marker genes
18. Save outputs

## QC rules

- Compute QC metrics before filtering
- Record thresholds
- Do not hardcode thresholds blindly
- Make mitochondrial detection configurable
- Save QC plots

## Normalization rules

- Preserve raw counts in adata.layers["counts"]
- Use config-driven parameters
- Avoid destructive overwrites

## Harmony rules

- Run after PCA, before neighbors
- Use only valid batch keys
- Never use biological condition as batch key
- Save embedding in adata.obsm["X_pca_harmony"]
- Use it in neighbors step
- Warn about overcorrection

## Clustering and markers

- Parameters must come from config
- Save marker tables
- Do not assign cell types without evidence

## Biology reviewer mode

- Use `prompts/biology_reviewer.md` when the user asks for a biology review of a run.
- Review only `biology_review_packet.md`, marker tables, QC summaries, and user-provided plots or screenshots unless the user explicitly allows deeper inspection.
- Do not access `data/` or inspect `.h5ad` files directly during biology review unless explicitly allowed.
- Provide tentative interpretations only, with marker evidence and uncertainty clearly separated.
- Recommend practical next checks such as continuing, changing Leiden resolution, changing the Harmony batch key, skipping Harmony, or inspecting specific clusters.
- Never recommend biological disease, diagnosis, treatment, phenotype, condition, or outcome variables as Harmony batch keys.

## Differential expression

- Distinguish marker vs condition DE
- Prefer replicate-aware methods
- Do not treat cells as independent samples

## Reproducibility

- Scripts over notebooks
- Config-driven pipeline
- Save parameters in adata.uns["analysis_params"]
- Save tables and figures

## Local validation

Allowed:

python -m py_compile scripts/*.py

Do NOT run full analysis locally.

## Server execution

Pipeline will run on server:

conda run -n scrna-agent python scripts/scrna_pipeline.py --config configs/pipeline.server.yaml

## Codex behavior

- Make minimal changes
- Do not fabricate data
- Do not run heavy computations
- Ask before risky operations

## Definition of done

Pipeline should produce:

- processed_with_harmony.h5ad
- QC summary CSV
- marker genes CSV
- UMAP plots
- saved parameters

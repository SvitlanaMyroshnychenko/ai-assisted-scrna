# scRNA-seq Codex Agent

Local template project for developing a reproducible single-cell RNA-seq analysis pipeline with Codex.

## Architecture

- Codex runs locally.
- Sensitive scRNA-seq data stays only on the remote server.
- Full analysis is executed manually on the remote server.
- Local files are templates for code, configs, notebooks, and SLURM scripts.

## Main files

- `AGENTS.md` — rules and instructions for Codex
- `configs/pipeline.template.yaml` — local template config
- `configs/samples.template.csv` — local template sample sheet
- `env/environment.yml` — conda environment
- `scripts/run_scanpy_pipeline.py` — main pipeline script, to be created
- `scripts/slurm_run_scanpy_pipeline.sh` — SLURM script, to be created
- `notebooks/01_view_results.ipynb` — result inspection notebook, to be created

## Safety

Do not commit:

- raw data
- processed `.h5ad`
- real server sample sheets
- private paths
- `.env`
- logs
- generated figures and tables
# Project Context

## Goal
Human-in-the-loop AI-assisted scRNA-seq/snRNA-seq workflow.

## Current Architecture
- deterministic Scanpy pipeline
- notebook orchestration workflow
- workflow_state.json state tracking
- notebook_helpers.py orchestration helpers
- staged review workflow

## Notebook Structure
0. Setup
1. Initialize Helpers
2. Initial Parameters
3. Write Config
4. QC
5. Preprocess
6. Harmony
7. Clustering
8. Markers
9. Review Workflow

## Important Design Decisions
- pipeline logic stays in scripts/scrna_pipeline.py
- notebook acts as workflow cockpit
- helpers live in scripts/notebook_helpers.py
- incremental refactors only
- no autonomous execution
- no LangChain/LangGraph rewrite
- stage-specific parameters are being moved closer to workflow stages

## Current Refactor Status
Completed:
- helper extraction
- workflow_state.json
- notebook cleanup
- write_config helper refactor
- QC parameter exposure

In Progress:
- stage-specific parameter UX
- notebook simplification

## Safety Constraints
- do not commit raw data
- do not commit .h5ad
- do not break pipeline behavior
- keep reviewer prompts unchanged unless explicitly requested
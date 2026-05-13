# Preprocess Best-Practice Checklist

Use this checklist when reviewing normalization, log transform, highly variable gene selection, scaling, and PCA.

## Goal

Decide whether preprocessing produced a technically reasonable representation for Harmony or clustering.

Do not assign cell types or interpret biological conditions at this stage.

## Expected Evidence

- Current normalization, HVG, scaling, and PCA parameters.
- Confirmation that raw counts are preserved.
- HVG summary and HVG plot.
- PCA scree plot.
- PCA scatter plots colored by QC metrics.
- Metadata-colored PCA plots when available.

## Review Questions

- Are raw counts preserved in `adata.layers["counts"]`?
- Is `adata.raw` present when expected?
- Is normalization/log transform conventional for this workflow?
- Does HVG selection look plausible across mean expression?
- Is the number of HVGs reasonable for the dataset complexity?
- If batches are heterogeneous, should batch-aware HVG selection be compared?
- Does PCA variance show a plausible elbow or long tail?
- Are PC1/PC2 driven by total counts, detected genes, mitochondrial percentage, or other technical metrics?
- Are enough PCs computed for downstream neighbor graph settings?

## Red Flags

- Raw counts not preserved.
- HVG selection dominated by technical artifacts.
- PCA strongly driven by QC metrics.
- Large unexplained PCA separation with no metadata diagnostics.
- Too few PCs computed for planned neighbors/clustering.
- Using all computed PCs downstream without justification when scree suggests fewer are appropriate.

## Acceptable Recommendations

- Proceed to Harmony.
- Skip Harmony and proceed to clustering if integration is not justified.
- Rerun preprocess with different HVG count, HVG flavor, batch-aware HVG selection, scaling, or PCA components.
- Add PCA plots colored by candidate batch/sample/QC metadata.

## Not Allowed

- Do not assign cell types.
- Do not infer disease effects.
- Do not request raw expression matrices.

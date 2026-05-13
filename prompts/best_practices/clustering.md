# Clustering Best-Practice Checklist

Use this checklist when reviewing neighbors, UMAP, and Leiden clustering.

## Goal

Decide whether the clustering is technically usable for marker analysis, or whether graph/clustering parameters should be revisited.

## Expected Evidence

- Neighbor parameters.
- Number of PCs used.
- UMAP random state.
- Leiden resolution.
- Number of clusters and cluster sizes.
- UMAP colored by Leiden.
- UMAP colored by QC metrics.
- UMAP colored by batch/sample/donor and known broad metadata when available.

## Review Questions

- Is the number of clusters plausible for the dataset size and debug/full context?
- Are there many tiny clusters?
- Are clusters driven by total counts, detected genes, mitochondrial percentage, or doublet score?
- Are clusters dominated by batch, donor, or sample?
- Does Harmony appear to preserve major structure while reducing technical segregation?
- Are neighbor graph parameters reasonable relative to PCA structure?
- Should Leiden resolution be increased, decreased, or kept?
- Can marker analysis proceed?

## Red Flags

- UMAP islands driven mostly by QC metrics.
- UMAP islands dominated by one batch/donor/sample.
- Many tiny clusters without coherent justification.
- Broad clusters likely mixing major cell classes.
- Clustering resolution produces fragmentation without marker evidence.
- Harmony appears to overcorrect or undercorrect.

## Acceptable Recommendations

- Proceed to marker analysis.
- Rerun clustering with lower or higher Leiden resolution.
- Adjust `neighbors_n_neighbors` or `neighbors_n_pcs`.
- Compare Harmony vs non-Harmony clustering.
- Inspect additional metadata overlays.

## Not Allowed

- Do not assign final cell types.
- Do not make disease or condition claims.
- Do not treat a debug subset clustering as final.

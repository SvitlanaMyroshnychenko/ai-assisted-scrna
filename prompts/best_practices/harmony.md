# Harmony Best-Practice Checklist

Use this checklist when reviewing Harmony or any PCA-level batch correction step.

## Goal

Decide whether Harmony was justified, whether the batch key is appropriate, and whether downstream clustering should use the adjusted representation.

## Expected Evidence

- Harmony enabled/disabled status.
- Batch key.
- Input basis and adjusted basis.
- Shape of `X_pca` and `X_pca_harmony`.
- Before/after diagnostic plots.
- Plots colored by batch key, QC metrics, and known broad annotations when available.

## Review Questions

- Is the batch key technical rather than biological or clinical?
- Does the batch key plausibly reflect donor, sequencing batch, chemistry, library prep, or other technical variation?
- Could the batch key also encode true biology?
- Does Harmony improve mixing without collapsing major biological structure?
- Are rare populations preserved?
- Do QC gradients still dominate after correction?
- Should clustering use `X_pca_harmony` or original `X_pca`?

## Red Flags

- Batch key is diagnosis, disease, condition, phenotype, treatment, cognition, clinical status, outcome, or similar.
- Major known cell classes collapse together after Harmony.
- Rare populations disappear or merge.
- Batch remains strongly segregated.
- Harmony is accepted without diagnostic plots.
- Donor-level correction is used without noting possible biological confounding.

## Acceptable Recommendations

- Proceed to clustering using `X_pca_harmony`.
- Proceed to clustering using original `X_pca`.
- Compare Harmony and non-Harmony clustering.
- Rerun Harmony with a different valid technical batch key.
- Skip Harmony if integration is not justified or overcorrection risk is high.

## Not Allowed

- Do not recommend biological, disease, condition, phenotype, diagnosis, treatment, cognition, clinical status, or outcome variables as batch keys.
- Do not assume Harmony improves results without evidence.
- Do not make disease interpretations from integration behavior.

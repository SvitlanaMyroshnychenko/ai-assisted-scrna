# Marker Best-Practice Checklist

Use this checklist when reviewing marker-gene results for clusters.

## Goal

Decide whether cluster marker programs are coherent enough for canonical marker review or manual annotation draft.

Marker analysis is evidence generation, not final annotation.

## Expected Evidence

- Marker summary.
- Marker table preview by cluster.
- Cluster sizes.
- Marker dotplot.
- Marker matrixplot.
- UMAP/clustering context when available.

## Review Questions

- Does each cluster have coherent top markers?
- Are markers specific, interpretable, and biologically plausible?
- Do related clusters share expected lineage markers?
- Are tiny clusters supported by coherent markers or likely unstable?
- Are any clusters dominated by stress, mitochondrial, ribosomal, hemoglobin, or dissociation signatures?
- Are mixed lineage markers present in the same cluster?
- Does marker evidence suggest Leiden resolution should change?
- Is canonical marker review/manual annotation justified?

## Red Flags

- Weak or nonspecific marker genes.
- Top markers dominated by technical/stress/RBC/ribosomal programs.
- Mixed unrelated lineage markers in one cluster.
- Tiny clusters with limited evidence.
- Marker programs inconsistent with UMAP or metadata context.
- Marker table interpreted as condition differential expression.

## Acceptable Recommendations

- Proceed to canonical marker review.
- Keep clusters exploratory until full run.
- Revisit tiny or mixed clusters.
- Rerun clustering at a different resolution.
- Add canonical marker panels or marker dotplots.

## Not Allowed

- Do not assign final labels automatically.
- Do not make disease or condition interpretations.
- Do not treat marker genes as replicate-aware condition DE.
- Do not request raw matrices.

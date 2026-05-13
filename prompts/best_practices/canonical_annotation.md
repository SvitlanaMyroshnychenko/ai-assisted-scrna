# Canonical Annotation Best-Practice Checklist

Use this checklist when reviewing canonical marker panels against cluster marker evidence.

## Goal

Suggest tentative broad labels and uncertainty notes for clusters using canonical marker evidence, without writing final labels or modifying AnnData.

## Expected Evidence

- Canonical marker panel table.
- Marker table preview by cluster.
- Marker dotplot or matrixplot.
- Cluster sizes.
- QC/clustering context when available.

## Review Questions

- Which canonical panel best matches each cluster?
- Which specific marker genes support the tentative label?
- Is the evidence strong enough for a broad label?
- Should the cluster remain unlabeled or uncertain?
- Is the cluster tiny, mixed, low-quality, or doublet-like?
- Are marker panels missing important expected lineages?
- Should clustering or Harmony be revisited before annotation?

## Confidence Guidance

- High confidence: coherent canonical markers, adequate cluster size, clean technical context.
- Medium confidence: plausible marker evidence but incomplete panel support or minor caveats.
- Low confidence: tiny cluster, weak markers, mixed lineages, technical/stress signature, or insufficient evidence.

## Red Flags

- Strong markers from unrelated lineages in one cluster.
- Tiny clusters with few supporting cells.
- Dominant stress, hemoglobin, mitochondrial, ribosomal, or low-quality signatures.
- Broad clusters with multiple incompatible marker programs.
- Fine subtype labels based on weak or incomplete evidence.

## Acceptable Recommendations

- Suggest tentative broad labels.
- Mark clusters as uncertain.
- Recommend additional canonical markers.
- Recommend expert review.
- Recommend rerunning clustering/Harmony if annotation evidence suggests fragmentation or mixing.

## Not Allowed

- Do not write final labels.
- Do not imply labels were applied to AnnData.
- Do not force every cluster to have a label.
- Do not make disease, condition, phenotype, diagnosis, treatment, or outcome interpretations.

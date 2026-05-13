# QC Best-Practice Checklist

Use this checklist when reviewing the QC stage.

## Goal

Decide whether the dataset is technically clean enough to proceed to preprocessing, or whether QC thresholds/diagnostics should be revisited.

Do not interpret cell types, disease biology, or condition effects at this stage.

## Expected Evidence

- QC filtering summary.
- Number and fraction of cells retained/removed.
- Current QC thresholds.
- Histograms or violin plots for detected genes, total counts, mitochondrial percentage, ribosomal percentage, and hemoglobin percentage when available.
- Scatter plot of counts vs detected genes colored by mitochondrial percentage.
- Doublet summary if doublet columns are available.
- Per-sample or per-batch QC summaries when available.

## Review Questions

- Are filtering thresholds explicitly recorded and biologically/technically justified?
- Does filtering appear too strict or too permissive?
- Are low-complexity cells removed without cutting into the main population?
- Is mitochondrial filtering appropriate for the tissue and chemistry?
- Does a high-count/high-gene tail remain after QC?
- Are high-count/high-gene cells enriched for doublet calls when doublet metadata exists?
- Is hemoglobin/RBC contamination visible?
- Is ribosomal percentage unusually high or clustered?
- Are removed cells balanced across samples/batches, or does one sample dominate the losses?

## Red Flags

- Very high fraction of cells removed without explanation.
- Very low fraction removed despite obvious low-quality tails.
- High mitochondrial tail retained after filtering.
- Strong high-count/high-gene tail with no doublet assessment.
- Zero doublets in a context where doublet detection should have been active.
- QC metrics strongly batch/sample-specific.
- Missing threshold rationale.

## Acceptable Recommendations

- Proceed to preprocessing.
- Proceed cautiously and inspect additional QC diagnostics.
- Rerun QC with adjusted thresholds.
- Add per-sample QC summaries before full run.
- Add threshold overlays or doublet/high-library-size diagnostics.

## Not Allowed

- Do not assign cell types.
- Do not make disease or condition interpretations.
- Do not request raw matrices.
- Do not recommend using biological/clinical metadata as correction variables.

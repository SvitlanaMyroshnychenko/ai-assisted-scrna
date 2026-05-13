Use prompts/biology_reviewer.md as your reviewer instructions.

Please review this single-cell run using only the packet content below and any plots/screenshots I provide. Do not access data/ or inspect h5ad files directly unless I explicitly allow it.


## Harmony Review Guidance

Please evaluate:

- Whether Harmony was justified for this dataset.
- Whether the selected batch key is appropriate.
- Whether the batch key might encode biological or clinical variation.
- Whether Harmony ran or was correctly skipped.
- Whether overcorrection risk should be checked before clustering.
- Whether clustering should use X_pca_harmony or original X_pca.

# Biology Review Packet

- Review step: `harmony`
- Run: `ec_debug_20260513_143103`
- Config: `configs/pipeline.server.yaml`

## Current Stage Parameters

### Harmony
- enabled: `True`
- batch_key: `individualID`
- basis: `X_pca`
- adjusted_basis: `X_pca_harmony`

## QC Summary
| sample_id   |   cells_before_filtering |   cells_after_filtering |   cells_removed |   min_genes |   min_cells |   max_pct_mt |   max_genes |   max_counts | mad_filter_enabled   |   log1p_total_counts_nmads |   log1p_n_genes_by_counts_nmads |   pct_counts_mt_nmads | mt_direction   |
|:------------|-------------------------:|------------------------:|----------------:|------------:|------------:|-------------:|------------:|-------------:|:---------------------|---------------------------:|--------------------------------:|----------------------:|:---------------|
| all         |                     5000 |                    4255 |             745 |         201 |          10 |            8 |         nan |          nan | True                 |                          5 |                               5 |                     3 | upper          |

## Harmony Summary
| harmony_enabled   | harmony_ran   | harmony_skipped   |   skip_reason | batch_key    |   n_batch_groups | batch_groups_preview                                                                                               | basis   | basis_shape   | adjusted_basis   | adjusted_basis_shape   |   risky_batch_key_warning |
|:------------------|:--------------|:------------------|--------------:|:-------------|-----------------:|:-------------------------------------------------------------------------------------------------------------------|:--------|:--------------|:-----------------|:-----------------------|--------------------------:|
| True              | True          | False             |           nan | individualID |               46 | R1924801; R8329066; R3292822; R6370138; R4260171; R9880904; R6312023; R4637403; R4174623; R5031238; ... (46 total) | X_pca   | (4255, 55)    | X_pca_harmony    | (4255, 55)             |                       nan |

## Latest AnnData Summary
- Latest h5ad: `harmony.h5ad`
- Cells: `4255`
- Genes: `15519`
- obs columns (showing first 25 of 87): `region`, `braaksc`, `inh.subtype`, `spanish`, `is.doublet`, `cts_mmse30_lv`, `individualIdSource`, `cogdx`, `individualID`, `cts_mmse30_first_ad_dx`, `hcluster`, `major.celltype`, `neuronal.layer`, `clinical_pathological_AD`, `lbl`, `U2`, `rind`, `sex`, `clinical_diagnosis`, `ceradsc`, `tspcol`, `col`, `hcelltype`, `full.exttype`, `neuronal.exttype`
- highly variable genes: `3050`
- raw present: `True`
- X_pca present: `True`
- X_pca shape: `(4255, 55)`
- X_pca_harmony present: `True`
- X_pca_harmony shape: `(4255, 55)`
- obsm keys: `X_pca`, `X_pca_harmony`
- layers: `counts`

## Marker Table
Marker table not available for this step.

## Figure Summary
Figure files not available for this step.

This packet intentionally excludes raw expression matrix data.

# Biology Review Packet

- Review step: `preprocess`
- Run: `ec_debug_20260513_143103`
- Config: `configs/pipeline.server.yaml`

## Current Stage Parameters

### Normalization
- target_sum: `10000`
- log1p: `True`

### Scaling
- max_value: `10`

### Highly Variable Genes
- n_top_genes: `3050`
- flavor: `seurat`
- batch_key: `None`

### PCA
- n_comps: `55`
- svd_solver: `arpack`
- random_state: `42`

### Debug
- max_cells: `5000`
- random_state: `42`

## QC Summary
| sample_id   |   cells_before_filtering |   cells_after_filtering |   cells_removed |   min_genes |   min_cells |   max_pct_mt |   max_genes |   max_counts | mad_filter_enabled   |   log1p_total_counts_nmads |   log1p_n_genes_by_counts_nmads |   pct_counts_mt_nmads | mt_direction   |
|:------------|-------------------------:|------------------------:|----------------:|------------:|------------:|-------------:|------------:|-------------:|:---------------------|---------------------------:|--------------------------------:|----------------------:|:---------------|
| all         |                     5000 |                    4255 |             745 |         201 |          10 |            8 |         nan |          nan | True                 |                          5 |                               5 |                     3 | upper          |

## Latest AnnData Summary
- Latest h5ad: `preprocessed.h5ad`
- Cells: `4255`
- Genes: `15519`
- obs columns (showing first 25 of 87): `region`, `braaksc`, `inh.subtype`, `spanish`, `is.doublet`, `cts_mmse30_lv`, `individualIdSource`, `cogdx`, `individualID`, `cts_mmse30_first_ad_dx`, `hcluster`, `major.celltype`, `neuronal.layer`, `clinical_pathological_AD`, `lbl`, `U2`, `rind`, `sex`, `clinical_diagnosis`, `ceradsc`, `tspcol`, `col`, `hcelltype`, `full.exttype`, `neuronal.exttype`
- highly variable genes: `3050`
- raw present: `True`
- X_pca present: `True`
- X_pca shape: `(4255, 55)`
- obsm keys: `X_pca`
- layers: `counts`

## Marker Table
Marker table not available for this step.

## Figure Summary
Review the attached contact sheet: `plots_contact_sheet_ec_debug_20260513_143103_preprocess.png`

This packet intentionally excludes raw expression matrix data.

# Biology Review Packet

- Review step: `qc`
- Run: `ec_debug_20260513_143103`
- Config: `configs/pipeline.server.yaml`

## Current Stage Parameters

### QC
- min_genes: `201`
- min_cells: `10`
- max_pct_mt: `8`
- max_genes: `None`
- max_counts: `None`
- mt_gene_prefix: `MT-`
- ribo_gene_prefixes: `['RPS', 'RPL']`
- hb_gene_pattern: `^HB[^(P)]`
- mad_filter.enabled: `True`
- mad_filter.log1p_total_counts_nmads: `5`
- mad_filter.log1p_n_genes_by_counts_nmads: `5`
- mad_filter.pct_counts_mt_nmads: `3`
- mad_filter.mt_direction: `upper`

### Debug
- max_cells: `5000`
- random_state: `42`

## QC Summary
| sample_id   |   cells_before_filtering |   cells_after_filtering |   cells_removed |   min_genes |   min_cells |   max_pct_mt |   max_genes |   max_counts | mad_filter_enabled   |   log1p_total_counts_nmads |   log1p_n_genes_by_counts_nmads |   pct_counts_mt_nmads | mt_direction   |
|:------------|-------------------------:|------------------------:|----------------:|------------:|------------:|-------------:|------------:|-------------:|:---------------------|---------------------------:|--------------------------------:|----------------------:|:---------------|
| all         |                     5000 |                    4255 |             745 |         201 |          10 |            8 |         nan |          nan | True                 |                          5 |                               5 |                     3 | upper          |

## Latest AnnData Summary
- Latest h5ad: `qc.h5ad`
- Cells: `4255`
- Genes: `15519`
- obs columns (showing first 25 of 87): `region`, `braaksc`, `inh.subtype`, `spanish`, `is.doublet`, `cts_mmse30_lv`, `individualIdSource`, `cogdx`, `individualID`, `cts_mmse30_first_ad_dx`, `hcluster`, `major.celltype`, `neuronal.layer`, `clinical_pathological_AD`, `lbl`, `U2`, `rind`, `sex`, `clinical_diagnosis`, `ceradsc`, `tspcol`, `col`, `hcelltype`, `full.exttype`, `neuronal.exttype`
- obsm keys: 
- layers: 
## Doublet Summary
- doublet call column: `is.doublet`
- called doublets: `0` / `4255` (`0.0%`)
- high-count/high-gene tail doublets (top 5% by either metric): `0` / `244` (`0.0%`)


## Marker Table
Marker table not available for this step.

## Figure Summary
Review the attached contact sheet: `plots_contact_sheet_ec_debug_20260513_143103_qc.png`

This packet intentionally excludes raw expression matrix data.

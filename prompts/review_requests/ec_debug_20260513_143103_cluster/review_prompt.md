Use prompts/biology_reviewer.md as your reviewer instructions.

Please review this single-cell run using only the packet content below and any plots/screenshots I provide. Do not access data/ or inspect h5ad files directly unless I explicitly allow it.

Please inspect the attached contact sheet image: `plots_contact_sheet_ec_debug_20260513_143103_cluster.png`.


## Clustering Review Guidance

Please evaluate:

- Whether UMAP structure looks plausible.
- Whether Leiden resolution appears too coarse or too granular.
- Whether clusters appear driven by QC metrics.
- Whether Harmony may have overcorrected or undercorrected.
- Whether neighbor graph parameters should be adjusted.
- Whether Leiden resolution should be increased or decreased.
- Whether marker analysis can safely proceed.

# Biology Review Packet

- Review step: `cluster`
- Run: `ec_debug_20260513_143103`
- Config: `configs/pipeline.server.yaml`

## Current Stage Parameters

### Neighbors
- n_neighbors: `13`
- n_pcs: `30`

### UMAP
- random_state: `42`

### Clustering
- method: `leiden`
- resolution: `0.5`
- key_added: `leiden`
- metadata_keys: `['leiden', 'individualID', 'major.celltype', 'hcelltype', 'region']`
- random_state: `42`

## QC Summary
| sample_id   |   cells_before_filtering |   cells_after_filtering |   cells_removed |   min_genes |   min_cells |   max_pct_mt |   max_genes |   max_counts | mad_filter_enabled   |   log1p_total_counts_nmads |   log1p_n_genes_by_counts_nmads |   pct_counts_mt_nmads | mt_direction   |
|:------------|-------------------------:|------------------------:|----------------:|------------:|------------:|-------------:|------------:|-------------:|:---------------------|---------------------------:|--------------------------------:|----------------------:|:---------------|
| all         |                     5000 |                    4255 |             745 |         201 |          10 |            8 |         nan |          nan | True                 |                          5 |                               5 |                     3 | upper          |

## Clustering Summary
| input_h5ad                                                                             |   cells |   genes | neighbors_use_rep   |   n_neighbors |   n_pcs |   umap_random_state | clustering_method   |   leiden_resolution | clustering_key   |   n_clusters | cluster_sizes_preview                                                                  | X_umap_present   | X_umap_shape   | X_pca_present   | X_pca_harmony_present   |
|:---------------------------------------------------------------------------------------|--------:|--------:|:--------------------|--------------:|--------:|--------------------:|:--------------------|--------------------:|:-----------------|-------------:|:---------------------------------------------------------------------------------------|:-----------------|:---------------|:----------------|:------------------------|
| /mnt/12tb_dsk3/svitlana/scrna-agent/results/ec_debug_20260513_143103/data/harmony.h5ad |    4255 |   15519 | X_pca_harmony       |            13 |      30 |                  42 | leiden              |                 0.5 | leiden           |           15 | 0:1109; 1:559; 10:102; 11:18; 12:17; 13:16; 14:13; 2:530; 3:514; 4:411; ... (15 total) | True             | (4255, 2)      | True            | True                    |

## Latest AnnData Summary
- Latest h5ad: `clustered.h5ad`
- Cells: `4255`
- Genes: `15519`
- obs columns (showing first 25 of 88): `region`, `braaksc`, `inh.subtype`, `spanish`, `is.doublet`, `cts_mmse30_lv`, `individualIdSource`, `cogdx`, `individualID`, `cts_mmse30_first_ad_dx`, `hcluster`, `major.celltype`, `neuronal.layer`, `clinical_pathological_AD`, `lbl`, `U2`, `rind`, `sex`, `clinical_diagnosis`, `ceradsc`, `tspcol`, `col`, `hcelltype`, `full.exttype`, `neuronal.exttype`
- highly variable genes: `3050`
- raw present: `True`
- X_pca present: `True`
- X_pca shape: `(4255, 55)`
- X_pca_harmony present: `True`
- X_pca_harmony shape: `(4255, 55)`
- X_umap present: `True`
- X_umap shape: `(4255, 2)`
- Leiden key present: `True`
- number of clusters: `15`
- cluster sizes preview: `0:1109; 1:559; 10:102; 11:18; 12:17; 13:16; 14:13; 2:530; 3:514; 4:411; ... (15 total)`
- obsm keys: `X_pca`, `X_pca_harmony`, `X_umap`
- layers: `counts`

## Marker Table
Marker table not available for this step.

## Figure Summary
Review the attached contact sheet: `plots_contact_sheet_ec_debug_20260513_143103_cluster.png`

This packet intentionally excludes raw expression matrix data.

# Biology Review Packet

- Review step: `markers`
- Run folder: `/mnt/12tb_dsk3/svitlana/scrna-agent/results/ec_debug_20260511_082408`
- Config path: `/mnt/12tb_dsk3/svitlana/scrna-agent/configs/pipeline.server.yaml`

## Current Parameters
- DATASET_ID: `EC`
- RUN_ID: `ec_debug`
- use_harmony: `True`
- harmony_batch_key: `individualID`
- neighbors_n_neighbors: `15`
- neighbors_n_pcs: `30`
- hvg_n_top_genes: `3000`
- hvg_flavor: `seurat`
- leiden_resolution: `0.5`
- random_state: `42`
- debug_max_cells: `5000`

## QC Summary
| sample_id   |   cells_before_filtering |   cells_after_filtering |   cells_removed |   min_genes |   min_cells |   max_pct_mt |   max_genes |   max_counts | mad_filter_enabled   |   log1p_total_counts_nmads |   log1p_n_genes_by_counts_nmads |   pct_counts_mt_nmads | mt_direction   |
|:------------|-------------------------:|------------------------:|----------------:|------------:|------------:|-------------:|------------:|-------------:|:---------------------|---------------------------:|--------------------------------:|----------------------:|:---------------|
| all         |                     5000 |                    4255 |             745 |         200 |          10 |            8 |         nan |          nan | True                 |                          5 |                               5 |                     3 | upper          |

## Marker Summary
| input_h5ad                                                                               |   cells |   genes | groupby   | method   |   n_groups | groups                                   |   n_genes_requested | used_raw   | marker_table_path                                                                               |
|:-----------------------------------------------------------------------------------------|--------:|--------:|:----------|:---------|-----------:|:-----------------------------------------|--------------------:|:-----------|:------------------------------------------------------------------------------------------------|
| /mnt/12tb_dsk3/svitlana/scrna-agent/results/ec_debug_20260511_082408/data/clustered.h5ad |    4255 |   15519 | leiden    | wilcoxon |         13 | 0; 1; 10; 11; 12; 2; 3; 4; 5; 6; 7; 8; 9 |                 100 | True       | /mnt/12tb_dsk3/svitlana/scrna-agent/results/ec_debug_20260511_082408/tables/cluster_markers.csv |

## Latest AnnData Summary
- Latest h5ad: `/mnt/12tb_dsk3/svitlana/scrna-agent/results/ec_debug_20260511_082408/data/markers.h5ad`
- Cells: `4255`
- Genes: `15519`
- obs columns (showing first 25 of 88): `region`, `braaksc`, `inh.subtype`, `spanish`, `is.doublet`, `cts_mmse30_lv`, `individualIdSource`, `cogdx`, `individualID`, `cts_mmse30_first_ad_dx`, `hcluster`, `major.celltype`, `neuronal.layer`, `clinical_pathological_AD`, `lbl`, `U2`, `rind`, `sex`, `clinical_diagnosis`, `ceradsc`, `tspcol`, `col`, `hcelltype`, `full.exttype`, `neuronal.exttype`
- highly variable genes: `3000`
- raw present: `True`
- X_pca present: `True`
- X_pca shape: `(4255, 50)`
- X_pca_harmony present: `True`
- X_pca_harmony shape: `(4255, 50)`
- X_umap present: `True`
- X_umap shape: `(4255, 2)`
- Leiden key present: `True`
- number of clusters: `13`
- cluster sizes preview: `0:1409; 1:560; 10:18; 11:17; 12:17; 2:530; 3:527; 4:411; 5:316; 6:218; ... (13 total)`
- obsm keys: `X_pca`, `X_pca_harmony`, `X_umap`
- layers: `counts`

## Marker Table
|   cluster | gene        |   rank |    score |   logfoldchange |         pval |     pval_adj |
|----------:|:------------|-------:|---------:|----------------:|-------------:|-------------:|
|         0 | ST18        |      1 | 51.5031  |         8.82482 | 0            | 0            |
|         0 | MBP         |      2 | 51.4136  |         4.91696 | 0            | 0            |
|         0 | CTNNA3      |      3 | 51.1581  |         6.63069 | 0            | 0            |
|         0 | PIP4K2A     |      4 | 50.901   |         4.52742 | 0            | 0            |
|         0 | SLC44A1     |      5 | 50.7912  |         4.513   | 0            | 0            |
|         1 | IQCJ-SCHIP1 |      1 | 36.1497  |         5.01948 | 3.76019e-286 | 5.83544e-282 |
|         1 | ST6GALNAC5  |      2 | 34.7078  |         4.85636 | 6.0127e-264  | 4.66556e-260 |
|         1 | DLGAP2      |      3 | 34.6267  |         4.68793 | 1.0019e-262  | 5.18283e-259 |
|         1 | KCNIP4      |      4 | 34.429   |         5.0912  | 9.28422e-260 | 3.60205e-256 |
|         1 | HS6ST3      |      5 | 34.0331  |         4.67498 | 7.217e-254   | 2.24001e-250 |
|         2 | DTNA        |      1 | 35.4606  |         4.30523 | 1.99471e-275 | 3.0956e-271  |
|         2 | SORBS1      |      2 | 34.9769  |         4.0697  | 5.04793e-268 | 3.91694e-264 |
|         2 | ADGRV1      |      3 | 34.9014  |         6.92192 | 7.07583e-267 | 3.66033e-263 |
|         2 | GPM6A       |      4 | 34.7937  |         4.14869 | 3.02382e-265 | 1.17317e-261 |
|         2 | PITPNC1     |      5 | 34.6933  |         4.93639 | 9.94191e-264 | 3.08577e-260 |
|         3 | DOCK8       |      1 | 35.6069  |         8.95486 | 1.09618e-277 | 5.67054e-274 |
|         3 | PLXDC2      |      2 | 35.4991  |         4.34759 | 5.07381e-276 | 1.57481e-272 |
|         3 | APBB1IP     |      3 | 34.127   |         9.35136 | 2.9344e-255  | 5.69238e-252 |
|         3 | ARHGAP24    |      4 | 33.988   |         6.24543 | 3.35505e-253 | 5.78522e-250 |
|         3 | SFMBT2      |      5 | 33.7633  |         4.64123 | 6.81252e-250 | 9.61123e-247 |
|         4 | TNR         |      1 | 33.347   |         5.94298 | 8.05267e-244 | 1.24969e-239 |
|         4 | DSCAM       |      2 | 33.2563  |         5.40077 | 1.65537e-242 | 1.28449e-238 |
|         4 | LHFPL3      |      3 | 33.0127  |         6.78518 | 5.34465e-239 | 2.76479e-235 |
|         4 | PCDH15      |      4 | 32.9242  |         6.0712  | 9.90089e-238 | 3.8413e-234  |
|         4 | PTPRZ1      |      5 | 32.7991  |         5.73927 | 6.06928e-236 | 1.88378e-232 |
|         5 | RALYL       |      1 | 26.7472  |         4.57731 | 1.3325e-157  | 1.26785e-153 |
|         5 | RYR2        |      2 | 26.7395  |         4.08005 | 1.63394e-157 | 1.26785e-153 |
|         5 | SH3GL2      |      3 | 25.5393  |         3.78966 | 7.21277e-144 | 3.73116e-140 |
|         5 | MCTP1       |      4 | 25.4179  |         3.79311 | 1.60108e-142 | 6.21178e-139 |
|         5 | FRMPD4      |      5 | 24.9201  |         3.82481 | 4.50442e-137 | 1.39808e-133 |
|         6 | ADARB2      |      1 | 24.3409  |         5.16756 | 7.24405e-131 | 1.1242e-126  |
|         6 | RGS12       |      2 | 23.0196  |         4.06777 | 2.96458e-117 | 2.30037e-113 |
|         6 | ERBB4       |      3 | 22.4505  |         3.87818 | 1.26629e-111 | 6.55051e-108 |
|         6 | SYNPR       |      4 | 21.3756  |         5.11939 | 2.2529e-101  | 8.7407e-98   |
|         6 | GRIK2       |      5 | 21.3569  |         3.77885 | 3.36082e-101 | 1.04313e-97  |
|         7 | GRIP1       |      1 | 17.7588  |         5.03182 | 1.47469e-70  | 2.28858e-66  |
|         7 | CNTNAP2     |      2 | 17.6469  |         4.36701 | 1.07581e-69  | 8.24711e-66  |
|         7 | PAM         |      3 | 17.6246  |         3.71765 | 1.59426e-69  | 8.24711e-66  |
|         7 | NXPH1       |      4 | 17.5612  |         5.60495 | 4.88619e-69  | 1.89572e-65  |
|         7 | KCNC2       |      5 | 17.4527  |         4.71647 | 3.28187e-68  | 9.54896e-65  |
|         8 | FGF13       |      1 | 15.7892  |         6.26415 | 3.68989e-56  | 5.72634e-52  |
|         8 | GRIP1       |      2 | 15.5306  |         5.6883  | 2.15463e-54  | 1.67189e-50  |
|         8 | FSTL5       |      3 | 15.4666  |         5.76375 | 5.82719e-54  | 3.01441e-50  |
|         8 | KAZN        |      4 | 15.4211  |         3.94453 | 1.18007e-53  | 4.57839e-50  |
|         8 | GRIN2A      |      5 | 15.1952  |         4.47097 | 3.80627e-52  | 1.18139e-48  |
|         9 | FLT1        |      1 |  7.11235 |         9.8087  | 1.14086e-12  | 1.7705e-08   |
|         9 | ABCB1       |      2 |  6.67297 |         8.29797 | 2.50684e-11  | 1.59787e-07  |
|         9 | MECOM       |      3 |  6.64227 |         7.16176 | 3.08887e-11  | 1.59787e-07  |
|         9 | ST6GALNAC3  |      4 |  6.34607 |         3.9633  | 2.20878e-10  | 6.08394e-07  |
|         9 | CLDN5       |      5 |  6.33934 |        12.2945  | 2.30756e-10  | 6.08394e-07  |
|        10 | UACA        |      1 |  6.79493 |         5.93333 | 1.08368e-11  | 1.68177e-07  |
|        10 | TPM1        |      2 |  6.58217 |         4.52151 | 4.63622e-11  | 2.9733e-07   |
|        10 | SLC6A20     |      3 |  6.50622 |        10.2167  | 7.70633e-11  | 2.9733e-07   |
|        10 | NTRK3       |      4 |  6.47652 |         4.31554 | 9.38638e-11  | 2.9733e-07   |
|        10 | NID1        |      5 |  6.45143 |         7.49218 | 1.10803e-10  | 2.9733e-07   |
|        11 | PDGFRB      |      1 |  7.01582 |         6.75437 | 2.28597e-12  | 1.90779e-08  |
|        11 | CALD1       |      2 |  7.00564 |         5.00196 | 2.45866e-12  | 1.90779e-08  |
|        11 | DLC1        |      3 |  6.45439 |         5.06884 | 1.08654e-10  | 5.62067e-07  |
|        11 | EPS8        |      4 |  6.39673 |         5.00503 | 1.58744e-10  | 6.15887e-07  |
|        11 | ATP1A2      |      5 |  6.3169  |         5.44212 | 2.66857e-10  | 7.44778e-07  |
|        12 | F13A1       |      1 |  7.10851 |        10.0262  | 1.17306e-12  | 1.82046e-08  |
|        12 | FMN1        |      2 |  6.69544 |         5.21045 | 2.15017e-11  | 1.66843e-07  |
|        12 | TBXAS1      |      3 |  6.31809 |         5.11879 | 2.64816e-10  | 1.07584e-06  |
|        12 | RBM47       |      4 |  6.31097 |         5.8612  | 2.77295e-10  | 1.07584e-06  |
|        12 | MRC1        |      5 |  6.17852 |         7.26009 | 6.4704e-10   | 2.00828e-06  |

## Figure Files
- `/mnt/12tb_dsk3/svitlana/scrna-agent/results/ec_debug_20260511_082408/figures/marker_dotplot_top_genes.png`

This packet intentionally excludes raw expression matrix data.
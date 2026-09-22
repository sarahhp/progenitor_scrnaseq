# Variation in beige and white adipose progenitors


# References
- data/other_papers/Human.MitoCarta3.0.csv
- data/gene_sets/msigdb.v2024.1.Hs.symbols.gmt


# ggplot themes
theme_classic(base_size=15) + 
  theme(axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth=0.25)))


# Colour schemes

## white vs beige
### primary: 
- colorRamp2(c(-2,0,2),c("#1F78B4","beige", "#FF7F00"))
- cyan: dichromat::colorschemes$BluetoDarkOrange.12
### discrete:
- cyan: dichromat::colorschemes$BluetoDarkOrange.12[c(10,2)]

### beige continuous:
scale_colour_distiller(palette="Oranges", direction=1)

## Subject colours
pal = c("#DA1819", "#00AED1", "#434D51","#EBB400","#691A93","#3AAD00")

## Timecourse
### White only
Day0 to day3
- pal = c("#B3B3B3","#85C2EA","#1F78B4")
- extra replicate :   "#B2DF8A"
 
### White to beige 

Day 3 white ->  day 0 -> day3 beige 

-  pal = c("#1F78B4", "#85C2EA","#B3B3B3", "#F18486","#E31A1C")
- extra replicate "#FDBF6F"

## Clusters

### Progenitor clusters (p0-7)
- prog_cols = c("#DA1819","#691A93","#EBB400" , "#434D51","#00AED1","#3AAD00","#FF7F00","#AC8AD0")
- TBC!

### unintegrated progenitor clusters (u0+)

- uclusters =  ggsci::pal_simpsons()(12)[c(1:7,9:12)]


### Beige versus white clusters b0-12 & t0-12
- bclusters = tvthemes::gravityFalls_pal()(13)


# Figure Index
Access to raw figures from the manuscript. Aesthetic formatting has been applied when mounting for the manuscript (e.g. figure size, axis titles) but no data has been changed.   

## Figure 1
### Thermogenic programs in beige adipocytes are subject-specific  

[Figure1B](analysis/bulk/figures/bulk_day15_GO/day15_GSEA_top2pfilt-1.pdf)
[Figure1C](analysis/bulk/figures/bulk_day15_UCP1_expr/UCP1.CIDEA.CITED1.PM20D1_RNA_expr_two_stats.pdf)
[Figure1D](analysis/bulk/figures/bulk_day15_mitopathways/MitoPathways_GSEA_top3filt-3.pdf)
[Figure1E](analysis/bulk/figures/bulk_day15_mitopathways/NonCanonicalTherm_heatmap.pdf)
[Figure1F](analysis/bulk/figures/bulk_day15_mitopathways/OXPHOS_complexes_heatmap.pdf)

## Figure 2
### Figure 2: Adipose progenitors display subject-specific transcriptomes by scRNA-seq. 

[Figure2B](figures/figure2_progenitors_initial/B_umap_subjects-1.pdf)
[Figure2C](figures/figure2_progenitors_initial/C_cluster_composition-1.pdf)
[Figure2D](figures/figure2_progenitors_initial/D_umap_clusters-1.pdf)
[Figure2E](figures/figure2_progenitors_initial/E_cluster_go_terms-1.pdf)
[Figure2F](figures/figure2_progenitors_initial/F_donor_GO_top3-1.pdf)

## Figure 3
### Integrated progenitor transcriptomes

[Figure3A](figures/figure3_progenitors_integrated/A_umap-1.pdf)
[Figure3B](figures/figure3_progenitors_integrated/marker_genes-1.pdf)
[Figure3C](figures/figure3_progenitors_integrated/cluster_go_unique-1.pdf)
[Figure3D](figures/figure3_progenitors_integrated/D_cluster_composition_boxplot-1.pdf)

## Figure 4

[Figure4B](figures/figure4_adipogenesis/umap-1.pdf)
[Figure4C](figures/figure4_adipogenesis/umap-2.pdf)
[Figure4D](figures/figure4_adipogenesis/unnamed-chunk-3-1.pdf)
[Figure4E](figures/figure4_adipogenesis/go_heatmap-1.pdf)
[Figure4E legend](figures/figure4_adipogenesis/heatmap_legend-1.pdf)

## Figure 5
### Beige versus white adipogenesis

[Figure5B](figures/figure5_beige_vs_white/B_wvb_proportion_test-1.pdf)
[Figure5C](figures/figure5_beige_vs_white/C_nomt_de_per_cluster-1.pdf)
[Figure5D](figures/figure5_beige_vs_white/D_umap_time-1.pdf)
[Figure5E](figures/figure5_beige_vs_white/E_wvb_proportion_test-1.pdf)
[Figure5F](figures/figure5_beige_vs_white/F_top_markers_include_all-1.pdf)
[Figure5H](figures/figure5_beige_vs_white/H_toptss_bvw-1.pdf)

## Figure 7
### ICAM1-associated progenitor and adipogenesis analyses

[Figure7A](figures/figure7_icam1/progenitors_rpca-1.pdf)
[Figure7B](figures/figure7_icam1/progenitors_rpca-2.pdf)
[Figure7C](figures/figure7_icam1/white_only_adipogenesis-1.pdf)
[Figure7D](figures/figure7_icam1/white_only_adipogenesis-2.pdf)

## Supplementary Figures

### Supplementary Figure S2
[Supplementary Figure S2](figures/supp_figS2_progenitor_integration_trial/all_plots-1.pdf)

### Supplementary Figures S3-S4
[Supplementary Figure S3](figures/supp_figS3_S4_progenitors/3A_clustree-1.pdf)
[Supplementary Figure S4](figures/supp_figS3_S4_progenitors/4C_main_phasescatter-1.pdf)

### Supplementary Figure S6
[Supplementary Figure S6A](figures/supp_figS6_adipogenesis/Apca-1.pdf)
[Supplementary Figure S6B](figures/supp_figS6_adipogenesis/Btop_markers_logfc_pct1-1.pdf)
[Supplementary Figure S6C](figures/supp_figS6_adipogenesis/Cemont_annot-1.pdf)

### Supplementary Figure S7
[Supplementary Figure S7B](figures/supp_figS7_beige_vs_white-gene_level/B_umap-1.pdf)
[Supplementary Figure S7E](figures/supp_figS7_beige_vs_white-gene_level/E_wvb_proportion_test-1.pdf)
[Supplementary Figure S7F](figures/supp_figS7_beige_vs_white-gene_level/F_upper_day1_volcano-1.pdf)

### Supplementary Figures S8-S9
[Supplementary Figure S8](figures/supp_figS8_S9_beige_vs_white-TSS_level/S8_A,B_pparg_isoforms-1.pdf)
[Supplementary Figure S9A](figures/supp_figS8_S9_beige_vs_white-TSS_level/S9_Amarker_gene_by_time-1.pdf)
[Supplementary Figure S9B](figures/supp_figS8_S9_beige_vs_white-TSS_level/S9_B_go_heatmap-1.pdf)
[Supplementary Figure S9C](figures/supp_figS8_S9_beige_vs_white-TSS_level/S9_C_cell_numbers_min10cells-1.pdf)
[Supplementary Figure S9D](figures/supp_figS8_S9_beige_vs_white-TSS_level/S9_D_upper_day1_volcano-1.pdf)
[Supplementary Figure S9E](figures/supp_figS8_S9_beige_vs_white-TSS_level/S9_E_isoform_types-1.pdf)

### Supplementary Figure S10
[Supplementary Figure S10](figures/supp_figS10_tss_integration_trial/all_plots-1.pdf)

### Supplementary Figure S11
[Supplementary Figure S11A](figures/supp_figS11_bvw_trajectory/A_integration_trial_clusters-2.pdf)
[Supplementary Figure S11B](figures/supp_figS11_bvw_trajectory/B_white_monocle_umap-1.pdf)
[Supplementary Figure S11C](figures/supp_figS11_bvw_trajectory/C_wvb_violin-1.pdf)


## Figure 6


[Figure 6E](figures/figure6_early_adipogenic_trajectories/E_RNAvelocity-1.pdf)
[Figure 6F velocity](figures/figure6_early_adipogenic_trajectories/Fbvw_prog_cols-1.pdf)
[Figure 6F monocle](figures/figure6_early_adipogenic_trajectories/Dmonocle_split-1.pdf)
[Figure 6G](figures/figure6_early_adipogenic_trajectories/Gvelo_violin-1.pdf)


# Analysis workflows

**Differentiated adipocytes:**

bulk RNA-seq (Figure 1)

a. [GO](bulk_day15_GO.html)
b. [DGE between subjects](bulk_day15_beige_DGE_between_subjects.html)
c. [UCP1 expression](bulk_day15_UCP1_expr.html)
d. [MitoPathways](bulk_day15_mitopathways.html)


**Adipocyte progenitor analysis:**

day 0 scRNA-seq (Figure 2 & 3)

1. [Initial QC and analysis (2)](progenitors_initial_seurat_analysis.html)
2. [GO terms (2e)](progenitors_initial_marker_gene_ontologies.html)
3. [Integration trial (S2a)](progenitors_integration_trial.html): Test 3 integration methods (RPCA, Harmony and FastMNN) and save the resulting latent representations and UMAP/TSNE diensionality reductions. 
4. [check other cluster numbers (RPCA; S2b)](progenitors/progenitors_rpca_clustree.html)
5. [RPCA integration (3)](progenitors_rpca_integrations.html) 
6. [RPCA clusters and markers (3b)](progenitors_rpca_clusters_and_markers.html)*
7. [RPCA cell cycle (Supp 3)](progenitors_rpca_cell_cycle.html)* <- check cluster ids


Other tests

* [sctranform test](progenitors_sctransform_test.html)
* [try regressing out UMI count; no effect](progenitors/initial_seurat_analysis_regress_out.html)
* [RPCA published annotation](progenitors_rpca_published_annotation.html)*
* [RPCA Emont annotation](progenitors_rpca_azimuth_emont_annotation.html)


**Early adipogenesis**

scRNA-seq 
day 0, 1 and 3 of differentiation
Figure 4

1. [White adipogenesis](white_only_downsample.html)
2. [Go terms](white_only_downsample_markers.html)
3. [Composition?]()


**Beige vs White **

scRNA-seq,
day 0, 1 and 3 of differentiation with vs without rosiglitazone 

Gene level (Figure 5B & C)

1. [Initial Analysis (Supp 5 A,B)](beige_vs_white_downsample.html)
2. [GO terms (Supp 5C,D)](beige_vs_white_downsample_markers_and_go.html)
3. [Beige vs white DE per time (Supp 5F)](beige_vs_white_downsample_DE_per_time.html)
4. [Cluster composition test (5B)](beige_vs_white_downsample_cluster_tests.html)
5. [Beige vs white DE per cluster (5C)](beige_vs_white_downsample_DE_per_cluster.html)


TSS-level (Figure 5D-I)

1. [Initial analysis (5D)](tss_bvw_initial.html)
2. [Cluster composition test (5E)](tss_bvw_initial_cluster_tests.html)
3. [GO terms (5F,G)](tss_bvw_initial_markers_and_go.html)
4. [DE per time (Supp 7D)](tss_bvw_initial_DE_per_time.html)
5. [DE per cluster (5H,I)](tss_bvw_initial_DE_per_cluster.html)

**Trajectory analysis **

TSS-level trajectory analysis
Figure 6 A-C

1. [Integration trial (Supp 8)](tss_integration_trial.html)
2. [FastMNN analysis](tss_bvw_fastmnn.html)
3. [Monocle analysis on fastmnn (6A)](tss_bvw_fastmnn_monocle_umap.html)
4. [Monocle + progenitor cells (6B-C)](tss_bvw_fastmnn_monocle_umap_plots.html)


Gene-level trajectory analysis
Figure 6 D-H

1. [Integration trial (6D & Supp S9A-C)](white_only_integration_trial.html)
2.[RNA velocity analysis on fastmnn (6E,F)](white_only_fastmnn_velociraptor_dynamical.html)
3. [RNA velocity plots (6F)](white_only_fastmnn_velociraptor_dynamical_plots.html)
4. [Monocle analysis on fastmnn (6F)](white_only_fastmnn_monocle_tsne.html)
5. [Progenitor markers (6G,H)](white_only_fastmnn_progenitor_trajectory_markers.html)


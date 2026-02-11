"""   """

__author__ = "Sarah Hazell Pickering (s.h.pickering@medisin.uio.no)"
__date__ = "2026-01-12"

RENDER_RMD = """MODULE="RStudio/2023.12.1+402-1-R-4.2.1";
if ! module -t list 2>&1 | grep -q "^${{MODULE}}$"; then
   echo "Loading module: $MODULE";
   module load "$MODULE"; fi
Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/usr/lib/rstudio/resources/app/bin/quarto/bin/tools'); 
            rmarkdown::render('{}')\" 
"""


rule main_figures:
    input:
        "figures/figure2_progenitors_initial.html",
        "figures/figure3_progenitors_integrated.html",
        "figures/figure4_adipogenesis.html",
        "figures/figure5_beige_vs_white.html",
        "figures/figure6_early_adipogenic_trajectories.html"


ODIR = "output/bulk"
ADIR = "analysis/bulk"

rule render_rmd:
    input:
        dge = ODIR + "/beige_DGE_between_subjects.tsv",
        rmd = ADIR + "/bulk_day15_RNAseq_expr.Rmd"
    output:
        report = ADIR + "/bulk_day15_RNAseq_expr.html"
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}" 

##-------------------------------------------------##
##    Mature adipocytes (day 15) from 6 subjects   ##
##                  bulk RNA-seq                   ##
##-------------------------------------------------##

# 
# rule bulk_DGE:
#     input:
#         storage.http("https://github.com/sarahhp/splicing_thermogenesis/raw/refs/heads/main/03limma/DGElist_ob_limma_filt.RData")
#     output:
#         expand("{dir}/beige_DGE_between_subjects.tsv",
#                 dir = ODIR),
#         expand("{dir}/bulk_day15_beige_DGE_between_subjects.html",
#                 dir=ADIR)
#     script:
#         ADIR + "/bulk_day15_beige_DGE_between_subjects.Rmd"

rule bulk_RNAseq_expr:
    input:
        expand("{dir}/beige_DGE_between_subjects.tsv",
            dir=ODIR)
    output:
        expand("{dir}/bulk_day15_RNAseq_expr.html",
            dir=ADIR)
    script:
        ADIR + "/bulk_day15_RNAseq_expr.Rmd"

##-------------------------------------------------##
##      Progenitors (day 0) from 6 subjects        ##
##                   sc RNA-seq                    ##
##-------------------------------------------------##

PRO_ADIR = "analysis/progenitors"
INDIR = "data/cellranger/count_no_introns/day0"
ODIR = "output/progenitors/initial"

rule progenitors_initial:
    ''' Progenitors scRNAseq initial seurat analysis
    Status: Complete.    
    '''
    input:
        in10x = expand("{dir}/{sample}/outs/filtered_feature_bc_matrix/{files}",
                        dir = INDIR,
                        sample = expand("day0_subjectS{s}", s=range(1,6)),
                        files = ["barcodes.tsv.gz","features.tsv.gz","matrix.mtx.gz"]),
        rmd = PRO_ADIR + "/progenitors_initial_seurat_analysis.Rmd"
    output:
        #bpcells = directory(expand("{dir}/bpcells/{sample}",
         #               dir = ODIR,
          #              sample =  expand("day0_subjectS{s}", s=range()),
        report = PRO_ADIR + "/progenitors_initial_seurat_analysis.html",
        #raw_rdata = ODIR + "/unfiltered_object.rds", #if first time running only
        rdata = ODIR + "/complete_analysis.rds",
        subset_rdata = ODIR + "/10%_complete_analysis.rds",
        marker_genes = ODIR + "/marker_genes.txt",
        donor_markers =  ODIR +"/donor_marker_genes.txt",
        cluster_info = ODIR + "/cluster_composition.txt",
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}"


rule pi_markers_and_go:
    '''Status: Knitted'''
    input:
        rdata = ODIR + "/complete_analysis.rds",
        marker_genes = ODIR + "/marker_genes.txt",
        donor_markers =  ODIR +"/donor_marker_genes.txt",
        gene_sets = expand("data/gene_sets/{file}",
                            file = ["msigdb.v2024.1.Hs.symbols.gmt",
                            "abbreviations.txt",
                            "capitalisations.txt"]),
        rmd = PRO_ADIR + "/progenitors_initial_markers_and_GO.Rmd"
    output:
        report = PRO_ADIR + "/progenitors_initial_markers_and_GO.html",
        go = ODIR + "/ORA_marker_genes.txt",
        donor_go = ODIR + "/ORA_donor_marker_genes.txt",
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}"

rule figure2:
    '''Status: Knitted'''
    input:
        rdata = ODIR + "/complete_analysis.rds",
        cluster_info = ODIR + "/cluster_composition.txt",
        go = ODIR + "/ORA_marker_genes.txt",
        donor_go = ODIR + "/ORA_donor_marker_genes.txt",
        rmd = "figures/figure2_progenitors_initial.Rmd"
    output:
        report = "figures/figure2_progenitors_initial.html",
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}"

## trial integration with different methods
INDIR = ODIR
TRIAL_DIR = "output/progenitors/integration_trial"

rule progenitors_integration_trial:
    ''' 
    Status: Complete
    Harmony, FastMNN and RPCA integrations
    '''
    input:
        rdata = INDIR + "/complete_analysis.rds",
        rmd = PRO_ADIR + "/progenitors_integration_trial.Rmd"
    output:
        report = PRO_ADIR + "/progenitors_integration_trial.html",
        rdata = TRIAL_DIR + "/complete_analysis.rds",
        subset_rdata = TRIAL_DIR + "/10%_complete_analysis.rds",
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}"
        
### scvii integration has an indepdent snakemake file
include: "run_scvi.py"


SCVI_DIR = "output/progenitors/scvi_integration"

rule progenitors_scvi:
    '''Status: Complete'''
    input:
        anndata = SCVI_DIR + "/complete_analysis.scviintegrated.h5ad",
        rmd = PRO_ADIR + "/progenitors_scvi_integration.Rmd"
    output:
        report = PRO_ADIR + "/progenitors_scvi_integration.html",
        rdata = SCVI_DIR + "/complete_analysis.rds",
        subset_rdata = SCVI_DIR + "/10%_complete_analysis.rds",
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}"  
        
rule supp_figS2_integration_trial:
    input:
        tri_integr = TRIAL_DIR + "/complete_analysis.rds",
        scvi = SCVI_DIR + "/complete_analysis.rds",
        rmd = "figures/supp_figS2_progenitor_integration_trial.Rmd"
    output:
        report = "figures/supp_figS2_progenitor_integration_trial.html"
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}"
    
       
#Trial clustering at different resolutions
rule progenitors_rpca_clustree:
    ''' 
    Status: Complete
    '''
    input:
        rdata = TRIAL_DIR + "/complete_analysis.rds",
        rmd = PRO_ADIR + "/progenitors_rpca_clustree.Rmd"
    output:
        report = PRO_ADIR + "/progenitors_rpca_clustree.html",
        rdata = TRIAL_DIR + "/multiple_resolutions.rds",
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}"
        
        
ODIR = "output/progenitors/rpca"

rule progenitors_rpca_integration:
    '''Starting from initial analysis apply optimal integration method 
    and clustering 
    Status: Complete
    '''
    input:
        rdata = INDIR + "/complete_analysis.rds",
        rmd = PRO_ADIR + "/progenitors_rpca_integration.Rmd"  
    output:
        report = PRO_ADIR + "/progenitors_rpca_integration.html",
        rdata = ODIR + "/complete_analysis.rds",
        subset_rdata = ODIR + "/10%_complete_analysis.rds",
        marker_genes = ODIR + "/marker_genes.txt",
        cluster_info = ODIR + "/cluster_composition.txt",
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}"
              
 
rule rpca_markers_and_go:
    '''Status: Complete'''
    input:
        rdata = ODIR + "/complete_analysis.rds",
        marker_genes = ODIR + "/marker_genes.txt",
        gene_sets = expand("data/gene_sets/{file}",
                            file = ["msigdb.v2024.1.Hs.symbols.gmt",
                            "abbreviations.txt",
                            "capitalisations.txt"]),
        rmd = PRO_ADIR + "/progenitors_rpca_markers_and_GO.Rmd"
    output:
        report = PRO_ADIR + "/progenitors_rpca_markers_and_GO.html",
        go = ODIR + "/ORA_marker_genes.txt",
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}"      

rule figure3:
    '''Status: Complete'''
    input:
        rdata = ODIR + "/complete_analysis.rds",
        marker_genes = ODIR + "/marker_genes.txt",
        go = ODIR + "/ORA_marker_genes.txt",
        cluster_info = ODIR + "/cluster_composition.txt",
        rmd = "figures/figure3_progenitors_integrated.Rmd"
    output:
        report = "figures/figure3_progenitors_integrated.html"
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}"  

rule supp_figS3_S4_progenitors:
    '''Status: Knitted'''
    input:
        mult = TRIAL_DIR + "/multiple_resolutions.rds",#clustree
        rdata = ODIR + "/complete_analysis.rds", #cell cycle
        cluster_info = ODIR + "/cluster_composition.txt",
        rmd = "figures/supp_figS3_S4_progenitors.Rmd"
        #unintergrated progenitor markers
    output:
        report = "figures/supp_figS3_S4_progenitors.html"
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}"  
        


##-------------------------------------------------##
##        Adipogenesis day 0, day 1 & day3         ##
##                   sc RNA-seq                    ##
##-------------------------------------------------##

ADIR = "analysis/adipogenesis"
INDIR = "data/cellranger/count_r2_only"
ODIR = "output/adipogenesis/initial_r2_only"

rule initial_r2_only:
    '''Status: complete'''
    input:
        in10x = expand("{dir}/{sample}/outs/filtered_feature_bc_matrix/{files}",
                        dir = INDIR,
                        sample = ["white_day3","white_day1_rep2","white_day1",
                                    "asc_day0",
                                    "beige_day1", "beige_day1_rep2","beige_day3"],
                        files = ["barcodes.tsv.gz","features.tsv.gz","matrix.mtx.gz"]),
        rmd = ADIR + "/adipogenesis_initial_r2_only.Rmd"
    output:
        bpcells = directory(expand("{dir}/bpcells/{sample}",
                        dir = ODIR,
                        sample = ["white_day3","white_day1_rep2","white_day1",
                                    "asc_day0",
                                    "beige_day1", "beige_day1_rep2","beige_day3"])),
        report = ADIR + "/adipogenesis_initial_r2_only.html",
        raw_rdata = ODIR + "/unfiltered_object.rds",
        rdata = ODIR + "/complete_analysis.rds",
        subset_rdata = ODIR + "/10%_complete_analysis.rds",
        marker_genes = ODIR + "/marker_genes.txt",
        GO_table =  ODIR +"/ORA_marker_genes.txt",
        cluster_info = ODIR + "/cluster_composition.txt",
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}"     

##-------------------------------------------------##
##      Adipogenesis - white only downsample       ##
##-------------------------------------------------##

INDIR = ODIR
ODIR = "output/adipogenesis/white_only_downsample"

rule white_only_downsample:
    '''Status: complete'''
    input:
        rdata = INDIR + "/complete_analysis.rds",
        rmd = ADIR + "/white_only_downsample.Rmd"
    output:
        report = ADIR + "/white_only_downsample.html",
        rdata = ODIR + "/complete_analysis.rds",
        subset_rdata = ODIR + "/10%_complete_analysis.rds",
        marker_genes = ODIR + "/marker_genes.txt",
        GO_table =  ODIR +"/ORA_marker_genes.txt",
        cluster_info = ODIR + "/cluster_composition.txt",
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}"  
        
#already migrated: clustree

rule white_only_markers:
    '''Status: complete'''
    input:
        rdata = ODIR + "/complete_analysis.rds",
        marker_genes = ODIR + "/marker_genes.txt",
        GO_table =  ODIR +"/ORA_marker_genes.txt",
        rmd = ADIR + "/white_only_downsample_markers_and_go.Rmd"
    output:
        report = ADIR + "/white_only_downsample_markers_and_go.html",
        heatmap = ODIR + "/ORA_marker_genes_expression_matrix.tsv"
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}"  
        
rule emont:
    '''status: Knitted'''
    input:
        rdata = ODIR + "/complete_analysis.rds",
        rmd = ADIR + "/white_only_downsample_emont.Rmd"
    output:
        annot = ODIR + "/emont/complete_analysis.emont.rds",
        report = ADIR + "/white_only_downsample_emont.html"
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}" 
        
# rule make_lazaref:
#     '''status: Migrated
#     Knitting not possible, just for graphing
#     '''
#     input:
#         anndata = "data/other_papers/lazarescu2025/lazarescu_subq_all.h5ad",
#         rmd = "analysis/publ_data/lazarescu2025/lazarescu_create_azimuth_ref_adipo_only.Rmd"
#     output:
#         idx = expand("{dir}/azimuth_reference/{file}",
#                         dir = "output/publ_data/lazarescu2025",
#                         idx = ["idx.annoy","ref.Rds"])
#     params:
#         cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
#     shell:
#         "{params.cmd}" 
        

# use rule render_rmd as lazarescu:
#      '''Status: Ran :)     Knitting not possible, just for graphing '''
#     input:
#         idx = expand("{dir}/azimuth_reference/{file}",
#                         dir = "output/publ_data/lazarescu2025",
#                         idx = ["idx.annoy","ref.Rds"]),
#         rmd = ADIR + "/white_only_downsample_downsample_lazarescu.Rmd"
#     output:
#         annot = ODIR + "/lazarescu/complete_analysis.annot.rds"
     
# rule miranda:
#     '''Status tbcreated'''
                
rule figure4:
    '''Status: Complete
    '''
    input:
        rdata = ODIR + "/complete_analysis.rds",#B&DUMAPs
        cluster_info = ODIR + "/cluster_composition.txt",#C
        heatmap = ODIR + "/ORA_marker_genes_expression_matrix.tsv",#Eheatmap
        lazar = ODIR + "/lazarescu/complete_analysis.annot.rds",
        annot = ODIR + "/emont/complete_analysis_emont.rds",
        rmd = "figures/figure4_adipogenesis.Rmd"
    output:
        report = "figures/figure4_adipogenesis.html",
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}"  
        
        
##-------------------------------------------------##
##       Beige vs white (in early adipogenesis)    ##
##-------------------------------------------------##


ADIR = "analysis/beige_vs_white"
ODIR = "output/beige_vs_white/downsample"

rule beige_vs_white_downsample:
    '''Status: Knitted'''
    input:
        rdata = INDIR + "/complete_analysis.rds",
        rmd = ADIR + "/beige_vs_white_downsample.Rmd"
    output:
        report = ADIR + "/beige_vs_white_downsample.html",
        rdata = ODIR + "/complete_analysis.rds",
        subset_rdata = ODIR + "/10%_complete_analysis.rds",
        marker_genes = ODIR + "/marker_genes.txt",
        GO_table =  ODIR +"/ORA_marker_genes.txt",
        cluster_info = ODIR + "/cluster_composition.txt",
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}"  

#rule bvw_markers_and_go:
# '''Status: Knitted; no ouput files'''

rule bvw_DE_per_time:
 '''Status: Complete'''
    input:
        rdata = ODIR + "/complete_analysis.rds",
        rmd = ADIR + "/beige_vs_white_downsample_DE_per_time.Rmd"
    output:
        sc_DE_per_time = expand("{dir}/DE_{time}_wvb.tsv",
            dir=ODIR,time=["day1","day3"]),
        report = ADIR + "/beige_vs_white_downsample_DE_per_time.html"
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}"  

rule bvw_cluster_tests:
    '''Status: Complete '''
    input:
        rdata = ODIR + "/complete_analysis.rds",
        cluster_info = ODIR + "/cluster_composition.txt",
        rmd = ADIR + "/beige_vs_white_downsample_cluster_tests.Rmd"
    output:
        report = ADIR + "/beige_vs_white_downsample_cluster_tests.html",
        result = ODIR + "/cluster_composition_test.tsv",
        cluster_info = ODIR + "/cluster_composition_condition_time.tsv"
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}"  

rule bvw_DE_per_cluster:
    '''Status: Complete '''
    input:
        rdata = ODIR + "/complete_analysis.rds",
        cluster_info = ODIR + "/cluster_composition_condition_time.tsv",
        rmd = ADIR + "/beige_vs_white_downsample_DE_per_cluster.Rmd"
    output:
        report = ADIR + "/beige_vs_white_downsample_DE_per_cluster.html",
        deg = ODIR + "/DE_per_cluster_white_vs_beige.tsv",
        go = ODIR + "/DE_per_cluster_white_vs_beige_GO_terms.tsv"
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}" 

rule supp_figure5:
    input:
        rdata = ODIR + "/complete_analysis.rds",#B UMAP
        marker_genes = ODIR + "/marker_genes.txt", #C
        GO_table =  ODIR +"/ORA_marker_genes.txt", #D
        cluster_info = ODIR + "/cluster_composition_test.tsv",#new E
        sc_DE_per_time = expand("{dir}/DE_{time}_wvb.tsv",
                                dir = ODIR, 
                                time = ["day1","day3"]), #F per time volcano
                                #G bulk per time
                                #H no. DEG per cluster
        rmd = "figures/supp_figS5_beige_vs_white-gene_level.Rmd"
    output:
        report = "figures/supp_figS5_beige_vs_white-gene_level.html"
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}"  
        
##-------------------------------------------------##
##      Beige vs White Integration (for fig6)      ##
##-------------------------------------------------##

INDIR = ODIR
TRIAL_DIR = "output/beige_vs_white/bvw_integration_trial"
        
rule bvw_integration_trial:
    '''Status:  Knitted'''
    input:
        rdata = INDIR + "/complete_analysis.rds",
        rmd = ADIR + "/bvw_integration_trial.Rmd"
    output:
        rdata = TRIAL_DIR + "/complete_analysis.rds",
        subset_rdata = TRIAL_DIR + "/10%_complete_analysis.rds",
        report = ADIR + "/bvw_integration_trial.html",
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}"  
        

SCVI_DIR = "output/beige_vs_white/bvw_scvi_integration"

rule bvw_scvi:
    '''Status: Created'''
    input:
        anndata = SCVI_DIR + "/complete_analysis.scviintegrated.h5ad",
        rmd = ADIR + "/bvw_scvi_integration.Rmd"
    output:
        report = ADIR + "/bvw_scvi_integration.html",
        rdata = SCVI_DIR + "/complete_analysis.rds",
        subset_rdata = SCVI_DIR + "/10%_complete_analysis.rds",
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}"  
        
rule supp_figS11_bvw_integration:
    input:
        tri_integr = TRIAL_DIR + "/complete_analysis.rds",
        scvi = SCVI_DIR + "/complete_analysis.rds",
        rmd = "figures/supp_figS11_bvw_integration_trial.Rmd"
    output:
        report = "figures/supp_figS11_bvw_integration_trial.html"
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}"
        
ODIR = "output/beige_vs_white/bvw_fastmnn"        
        
rule bvw_fastmnn:
    '''Status:  created'''
    input:
        rdata = INDIR + "/complete_analysis.rds",
        rmd = ADIR + "/bvw_fastmnn_integration.Rmd"
    output:
        rdata = ODIR + "/complete_analysis.rds",
        marker_genes = ODIR + "/marker_genes.txt",
        GO_table =  ODIR +"/ORA_marker_genes.txt",
        cluster_info = ODIR + "/cluster_composition.tsv",
        report = ADIR + "/bvw_fastmnn_integration.html",
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}" 
    
## Decide to use UMAp or tsNE
rule bvw_monocle:
    '''Status:  Knitted'''
    input:
        rdata = ODIR + "/complete_analysis.rds",
        rmd = ADIR + "/bvw_fastmnn_monocle_umap.Rmd"
    output:
        white = ODIR + "/white_monocle_complete_analysis_umap.rds",
        beige = ODIR + "/beige_monocle_complete_analysis_umap.rds",
        report = ADIR + "/bvw_fastmnn_monocle_umap.html",
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}" 

HARMONY ="output/beige_vs_white/bvw_harmony" 
use rule bvw_monocle as bvw_monocle_harmony with:
    input:
        rdata = TRIAL_DIR + "/complete_analysis.rds",
        rmd = ADIR + "/bvw_harmony_monocle_umap.Rmd"
    output:
        white = HARMONY + "/white_monocle_complete_analysis_umap.rds",
        beige = HARMONY + "/beige_monocle_complete_analysis_umap.rds",
        report = ADIR + "/bvw_harmony_monocle_umap.html",
        
##-------------------------------------------------##
##                  TSS-level analysis             ##
##              of beige vs white data             ##
##-------------------------------------------------##

INDIR = "data/scafe/count/day0_to_day3_rep1"
TSS_DIR = "output/beige_vs_white/tss_bvw_initial"

rule tss_bvw_initial:
    '''Status: Complete'''
    input:
        raw_counts = expand("{dir}/{sample}/matrix/{files}",
                        dir = INDIR,
                        sample = ["white_day3","white_day1",
                                    "asc_day0",
                                    "beige_day1","beige_day3"],
                        files = ["barcodes.tsv","genes.tsv","matrix.mtx"]),
        rmd = ADIR + "/tss_bvw_initial.Rmd"
    output:
        #raw_rdata = TSS_DIR + "/unfiltered_object.rds", #if first time running only
        report = ADIR + "/tss_bvw_initial.html",
        rdata = TSS_DIR + "/complete_analysis.rds",
        subset_rdata = TSS_DIR + "/10%_complete_analysis.rds",
        marker_genes = TSS_DIR + "/marker_genes.txt",
        GO_table =  TSS_DIR +"/ORA_marker_genes.txt",
        cluster_info = TSS_DIR + "/cluster_composition.tsv",
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}"  
        
# rule tss_markers_and_go:
#     '''Status: Knitting; no output files but good reference for PPARG tss plots and GOheatmap'''
#  
# rule tss_DE_per_time:
#     '''Status: Migrated; not in figures but background info'''
        
rule tss_cluster_tests:
    '''Status: Complete'''
    input:
        rdata = TSS_DIR + "/complete_analysis.rds",
        cluster_info = TSS_DIR + "/cluster_composition.tsv",
        rmd = ADIR + "/tss_bvw_initial_cluster_tests.Rmd",
    output:
        report = ADIR + "/tss_bvw_initial_cluster_tests.html",
        cluster_info = TSS_DIR + "/cluster_composition_condition_time.tsv",
        result = TSS_DIR + "/cluster_composition_test.tsv",
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}" 

rule tss_DE_per_cluster:
    '''Status: Complete '''
    input:
        rdata = TSS_DIR + "/complete_analysis.rds",
        cluster_info = TSS_DIR + "/cluster_composition_condition_time.tsv",
        rmd = ADIR + "/tss_bvw_initial_DE_per_cluster.Rmd"
    output:
        report = ADIR + "/tss_bvw_initial_DE_per_cluster.html",
        deg = TSS_DIR + "/DE_per_cluster_wvb.tsv",
        go = TSS_DIR + "/DE_per_cluster_wvb_GO_terms.tsv"
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}" 

ODIR = "output/beige_vs_white/downsample"

rule figure5:
    input:
        #Beige vs white (gene level)
        bvw_result = ODIR + "/cluster_composition_test.tsv",#B cluster comp bvw
        bvw_deg = ODIR + "/DE_per_cluster_white_vs_beige.tsv",#C wvb genes per cluster
        #Tss level (beige vs white)
        rdata = TSS_DIR + "/complete_analysis.rds",#D
        result = TSS_DIR + "/cluster_composition_test.tsv",#E
        marker_genes = TSS_DIR + "/marker_genes.txt",#F 
        cluster_info = TSS_DIR + "/cluster_composition_condition_time.tsv", #for cluster annot
        GO_table =  TSS_DIR +"/ORA_marker_genes.txt",#G 
        tdeg = TSS_DIR + "/DE_per_cluster_wvb.tsv",#H
        go = TSS_DIR + "/DE_per_cluster_wvb_GO_terms.tsv",#I
        rmd = "figures/figure5_beige_vs_white.Rmd"
    output:
        report = "figures/figure5_beige_vs_white.html"
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}" 

##-------------------------------------------------##
##             Progenitor trajectories             ##
##              in early adipogenesis              ##
##-------------------------------------------------##

ODIR = "output/beige_vs_white/tss_integration_trial"

rule tss_integration_trial:
    '''Status: Complete
    Supp FigS8'''
    input:
        rdata = TSS_DIR + "/complete_analysis.rds",
        rmd = ADIR + "/tss_integration_trial.Rmd"
    output:
        rdata = ODIR + "/complete_analysis.rds",
        subset_rdata = ODIR + "/10%_complete_analysis.rds",
        report = ADIR + "/tss_integration_trial.html",
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}"
        
        
rule supp_figS10:
    '''tss_integration_trial'''
    input:
        rdata = ODIR + "/complete_analysis.rds",
        rmd = "figures/supp_figS10_tss_integration_trial.Rmd"
    output:
        report = "figures/supp_figS10_tss_integration_trial.html"
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}"
        
        
ODIR = "output/beige_vs_white/tss_fastmnn"
        
rule tss_fastmnn_integration:
    '''Status:  Knitting...'''
    input:
        rdata = TSS_DIR + "/complete_analysis.rds",
        rmd = ADIR + "/tss_fastmnn_integration.Rmd"
    output:
        rdata = ODIR + "/complete_analysis.rds",
        marker_genes = ODIR + "/marker_genes.txt",
        GO_table =  ODIR +"/ORA_marker_genes.txt",
        cluster_info = ODIR + "/cluster_composition.tsv",
        report = ADIR + "/tss_fastmnn_integration.html",
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}" 
        
rule tss_fastmnn_monocle:
    '''Status:  Complete'''
    input:
        rdata = ODIR + "/complete_analysis.rds",
        rmd = ADIR + "/tss_fastmnn_monocle_umap.Rmd"
    output:
        white = ODIR + "/white_monocle_complete_analysis_umap.rds",
        beige = ODIR + "/beige_monocle_complete_analysis_umap.rds",
        report = ADIR + "/tss_fastmnn_monocle_umap.html",
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}" 
        
        
rule tss_monocle_plots:
    input:
        white = ODIR + "/white_monocle_complete_analysis_umap.rds",
        beige = ODIR + "/beige_monocle_complete_analysis_umap.rds",
        rdata = ODIR + "/complete_analysis.rds",
        rmd = ADIR + "/tss_fastmnn_monocle_umap_plots.Rmd"
    output: 
        pseudotime = ODIR + "/monocle_pseudotime.tsv",
        seurat_white = ODIR + "/white_monocle_pseudotime_seurat.rds",
        seurat_beige = ODIR + "/beige_monocle_pseudotime_seurat.rds",
        report = ADIR + "/tss_fastmnn_monocle_umap_plots.html",
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}" 
   
BVW_DIR = "output/beige_vs_white/bvw_fastmnn"     
#Snakemake 6.0 and later, it is possible to inherit from previously defined rules
use rule tss_monocle_plots as bvw_monocle_plots with:
    input:
        white = BVW_DIR + "/white_monocle_complete_analysis_umap.rds",
        beige = BVW_DIR + "/beige_monocle_complete_analysis_umap.rds",
        rdata = BVW_DIR + "/complete_analysis.rds",
        rmd = ADIR + "/bvw_fastmnn_monocle_umap_plots.Rmd"
    output: 
        pseudotime = BVW_DIR + "/monocle_pseudotime.tsv",
        seurat_white = BVW_DIR + "/white_monocle_pseudotime_seurat.rds",
        seurat_beige = BVW_DIR + "/beige_monocle_pseudotime_seurat.rds",
        report = ADIR + "/bvw_fastmnn_monocle_umap_plots.html",
        
#include: "run_g2g.py"    

rule figure6:
    input:
        beige = ODIR + "/beige_monocle_complete_analysis_umap.rds", #A&B Trajectory TSSs
        pseudotime = ODIR + "/monocle_pseudotime.tsv", #C Violins
        rmd = "figures/figure6_early_adipogenic_trajectories.Rmd",
        gwhite = BVW_DIR + "/white_monocle_complete_analysis_umap.rds", #F? or D
        gbeige = BVW_DIR + "/beige_monocle_complete_analysis_umap.rds",
        gpseudotime = BVW_DIR + "/monocle_pseudotime.tsv",
    output:
        report = "figures/figure6_early_adipogenic_trajectories.html"
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}"        

SUPP_FIGS = [
     #"figures/supp_figS1_progenitors_unintegrated.html",
     "figures/supp_figS2_progenitor_integration_trial.html",
     "figures/supp_figS3_S4_progenitors.html",
    # "figures/supp_figS6_beige_vs_white-gene_level.html",
     
     "figures/supp_figS10_tss_integration_trial.html",
     #"figures/supp_figS11_bvw_integration_trial.html",
     
]
rule supp_figure:
     input:
         SUPP_FIGS
         
rule all:
    input: 
        "figures/figure2_progenitors_initial.html",
        "figures/figure3_progenitors_integrated.html",
        "figures/figure4_adipogenesis.html",
        "figures/figure5_beige_vs_white.html",
        "figures/figure6_early_adipogenic_trajectories.html",
        SUPP_FIGS









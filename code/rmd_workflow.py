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
 #shell: module load RStudio/2023.12.1+402-1-R-4.2.1
 #or check whther loaded: RStudio/2023.12.1+402-1-R-4.2.1 (L)

##-------------------------------------------------##
##    Mature adipocytes (day 15) from 6 subjects   ##
##                  bulk RNA-seq                   ##
##-------------------------------------------------##

ODIR = "output/bulk"
ADIR = "analysis/bulk"
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
    input:
        in10x = expand("{dir}/{sample}/outs/filtered_feature_bc_matrix/{files}",
                        dir = INDIR,
                        sample = expand("day0_subjectS{s}", s=range(1,6)),
                        files = ["barcodes.tsv.gz","features.tsv.gz","matrix.mtx.gz"]),
        rmd = PRO_ADIR + "/progenitors_initial_seurat_analysis.Rmd"
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    output:
        #bpcells = directory(expand("{dir}/bpcells/{sample}",
         #               dir = ODIR,
          #              sample =  expand("day0_subjectS{s}", s=range()),
        report = PRO_ADIR + "/progenitors_initial_seurat_analysis.html",
        #raw_rdata = ODIR + "unfiltered_object.rds", #if first time running only
        rdata = ODIR + "/complete_analysis.rds",
        subset_rdata = ODIR + "/10%_complete_analysis.rds",
        marker_genes = ODIR + "/marker_genes.txt",
        donor_markers =  ODIR +"/donor_marker_genes.txt",
        cluster_info = ODIR + "/cluster_composition.txt",
    shell:
        "{params.cmd}"


rule pi_markers_and_go:
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

## trial integration with different methods
INDIR = ODIR
ODIR = "output/progenitors/integration_trial"

rule progenitors_integration_trial:
    ''' 
    Status: knitted
    Harmony, FastMNN and RPCA integrations
    '''
    input:
        rdata = INDIR + "/complete_analysis.rds",
        rmd = PRO_ADIR + "/progenitors_integration_trial.Rmd"
    output:
        report = PRO_ADIR + "/progenitors_integration_trial.html",
        rdata = ODIR + "/complete_analysis.rds",
        subset_rdata = ODIR + "/10%_complete_analysis.rds",
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}"
        
### scvii integration has an indepdent snakemake file
include: "run_scvi.py"
#check the way I did it last time

SCVI_DIR = "output/progenitors/scvi_integration"

rule progenitors_scvi:
    input:
        anndata = SCVI_DIR + "/complete_analysis.scviintegrated.h5ad",
        rmd = PRO_ADIR + "/progenitors_scvi_integration.Rmd"
    output:
        report = PRO_ADIR + "/progenitors_integration_trial.html",
        rdata = SCVI_DIR + "/complete_analysis.rds",
        subset_rdata = SCVI_DIR + "/10%_complete_analysis.rds",
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}"  
       
#Trial clustering at different resolutions
rule progenitors_rpca_clustree:
    ''' 
    Status: knitted
    Harmony, FastMNN and RPCA integrations
    '''
    input:
        rdata = ODIR + "/complete_analysis.rds",
        rmd = PRO_ADIR + "/progenitors_rpca_clustree.Rmd"
    output:
        report = PRO_ADIR + "/progenitors_rpca_clustree.html",
        rdata = ODIR + "/multiple_resolutions.rds",
    params:
        cmd = lambda wildcards, input: RENDER_RMD.format(input.rmd)
    shell:
        "{params.cmd}"
        
        
ODIR = "output/progenitors/rpca"

rule progenitors_rpca_integration:
    '''Starting from initial analysis apply optimal integration method 
    and clustering 
    Status: knitted
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
              
        

##-------------------------------------------------##
##        Adipogenesis day 0, day 1 & day3         ##
##                   sc RNA-seq                    ##
##-------------------------------------------------##

ADIR = "analysis/adipogenesis"
INDIR = "data/cellranger/count_r2_only"
ODIR = "output/adipogenesis/initial_r2_only"

rule initial_r2_only:
    input:
        in10x = expand("{dir}/{sample}/outs/filtered_feature_bc_matrix/{files}",
                        dir = INDIR,
                        sample = ["white_day3","white_day1_rep2","white_day1",
                                    "asc_day0",
                                    "beige_day1", "beige_day1_rep2","beige_day3"],
                        files = ["barcodes.tsv.gz","features.tsv.gz","matrix.mtx.gz"]),
    output:
        bpcells = directory(expand("{dir}/bpcells/{sample}",
                        dir = ODIR,
                        sample = ["white_day3","white_day1_rep2","white_day1",
                                    "asc_day0",
                                    "beige_day1", "beige_day1_rep2","beige_day3"])),
        report = ADIR + "/adipogenesis_initial_r2_only.html",
        raw_rdata = ODIR + "unfiltered_object.rds",
        rdata = ODIR + "/complete_analysis.rds",
        subset_rdata = ODIR + "/10%_complete_analysis.rds",
        marker_genes = ODIR + "/marker_genes.txt",
        GO_table =  ODIR +"/ORA_marker_genes.txt",
        cluster_info = ODIR + "/cluster_composition.txt",
    script:
        ADIR + "/adipogenesis_initial_r2_only.Rmd"

INDIR = ODIR
ODIR = "output/adipogenesis/white_only_initial"
     
rule white_only_initial:
    ''' Not yet migrated '''
    input:
        rdata = INDIR + "/complete_analysis.rds"
    output:
        report = ADIR + "/white_only_initial.html",
        raw_rdata = ODIR + "unfiltered_object.rds",
        rdata = ODIR + "/complete_analysis.rds",
        subset_rdata = ODIR + "/10%_complete_analysis.rds",
        marker_genes = ODIR + "/marker_genes.txt",
        GO_table =  ODIR +"/ORA_marker_genes.txt",
        cluster_info = ODIR + "/cluster_composition.txt",
    script:
        ADIR + "/white_only_initial.Rmd"

##-------------------------------------------------##
##      Adipogenesis - white only downsample       ##
##-------------------------------------------------##

#INDIR = ODIR
ODIR = "output/adipogenesis/white_only_downsample"

rule white_only_downsample:
    input:
        rdata = INDIR + "/complete_analysis.rds"
    output:
        report = ADIR + "/white_only_downsample.html",
        raw_rdata = ODIR + "unfiltered_object.rds",
        rdata = ODIR + "/complete_analysis.rds",
        subset_rdata = ODIR + "/10%_complete_analysis.rds",
        marker_genes = ODIR + "/marker_genes.txt",
        GO_table =  ODIR +"/ORA_marker_genes.txt",
        cluster_info = ODIR + "/cluster_composition.txt",
    script:
        ADIR + "/white_only_downsample.Rmd"
        
#already migrated: clustree           
                
                
                


ODIR = "../output/bulk"
ADIR = "../analysis/bulk"
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
##        Adipogenesis day 0, day 1 & day3         ##
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
        
rue             
                
                
                

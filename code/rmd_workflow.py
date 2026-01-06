
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

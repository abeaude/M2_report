source("utils.R")
sourceFolder("Functionnal Enrichment")

library(future)
plan(multicore)


  # To avoid packages output startup useless
suppressPackageStartupMessages(library(get_db('mouse')['db'],character.only = TRUE))
path <- "/home/aurel/test/"

load("/home/aurel/Nextcloud/Documents/Cours/M2_GENIOMHE/Stage/Project/NF1-MPNST/NF1-MPNST_20190212/ANALYSIS/DEFAULT/NF1-MPNST_20190212_SCT_DEFAULT_pca100.rda")
run_multiple_DE(sobj = sobj, path = path)

run_DE(sobj = sobj, clusters = NULL, DE_type = "1vsAll",testToUse = NULL, method = c('Limma_Voom','Limma_Trend','EdgeR_QL'), first_group = NULL, second_group = NULL, path = "/home/aurel/Nextcloud/Documents/Cours/M2_GENIOMHE/Stage/Project/NF1-MPNST/NF1-MPNST_20190212/", batch = NULL,.parallel = FALSE)

# or
run_mutliple_DE(path = path, sobj = sobj)

stat_DE(path, species = "mouse")
run_functionnal_enrichment(path, pval = 0.05, thresholdLogFC = 0.58, thresholdpct = 5, species = "mouse", category_msigdb = c("H","C6"), method = c("WikiPathway", "msigdb", "GOterms", "ReactomePA"), gmt_date = NULL, gmt_folder = "/home/aurel/Nextcloud/Documents/Cours/M2_GENIOMHE/Stage/RESSOURCES/gmt_file", enrichment = TRUE, GSEA = TRUE, comparison = TRUE)

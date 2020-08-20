## Interaction
source("utils.R")
sourceFolder("Interaction/")

load("/home/aurel/Nextcloud/Documents/Cours/M2_GENIOMHE/Stage/Project/NF1-MPNST/NF1-MPNST_20190212/ANALYSIS/DEFAULT/NF1-MPNST_20190212_SCT_DEFAULT_pca100.rda")

path <- "/home/aurel/Nextcloud/Documents/Cours/M2_GENIOMHE/Stage/Project/NF1-MPNST/NF1-MPNST_20190212"
fs::path(path,"INTERACTION") %>% fs::dir_create()

CCInx <- run_CCInx(path,test = NULL, clusters = levels(Seurat::Idents(sobj)), species = "mouse")
saveRDS(CCInx, fs::path(path,"INTERACTION","CCInx", "results",ext = "rds"))
library(shiny)
CCInx::ViewCCInx(CCInx)

cellphoneDB <- run_cellphonedb(sobj, "gene_name", "mouse", path = path, filename = "test_phoneDB")
        
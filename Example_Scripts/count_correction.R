source("utils.R")
source("SoupX.R")

# Load raw matrix from cell ranger
raw_counts <- Seurat::Read10X("/media/aurel/ac3a5c6f-745e-49f7-9711-17b88dccbbe2/Nextcloud/Documents/Cours/M2_GENIOMHE/Stage/Project/NF1-MPNST/NF1-MPNST_20190212/counts/raw_feature_bc_matrix/")
# Load FILTERED
load("/media/aurel/ac3a5c6f-745e-49f7-9711-17b88dccbbe2/Nextcloud/Documents/Cours/M2_GENIOMHE/Stage/Project/NF1-MPNST/NF1-MPNST_20190212/ANALYSIS/DEFAULT/NF1-MPNST_20190212_SCT_DEFAULT_pca100.rda")


path <- '/media/aurel/ac3a5c6f-745e-49f7-9711-17b88dccbbe2/Nextcloud/Documents/Cours/M2_GENIOMHE/Stage/Project/NF1-MPNST/NF1-MPNST_20190212/'

# plot are stored in subfolder SOUPX in path
out <- remove_contamination(sobj, raw_counts,path)

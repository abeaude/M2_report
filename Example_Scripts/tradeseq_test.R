library(tradeSeq)
library(slingshot)

sds <- readRDS("/home/aurel/test/TRAJECTORIES/idents_1-4-7-8-12/slingshot/UMAP/trajectory_result_slingshot_dataset.rds")
sce <- readRDS("/home/aurel/test/TRAJECTORIES/idents_1-4-7-8-12/slingshot/data_for_trajectories.rds")
counts <- SingleCellExperiment::counts(sce) %>% as.matrix()

# find the number of knots
find_knots <- evaluateK(counts = counts, sds = sds)

# Used for paralleilzation
BPPARAM <- BiocParallel::bpparam()
# specify number of core to use
BPPARAM$workers <- 4

# running the fitGAM function
sce <- fitGAM(counts = counts, sds = sds,
            nknots = 6, verbose = TRUE,  sce = TRUE, parallel = FALSE)

assoRes <- associationTest(sce)
startRes <- startVsEndTest(sce)
oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce)[oStart[1]]

# Need multiple lineage
if(ncol(slingCurveWeights(sds)) > 1) {
  endRes <- diffEndTest(sce)
  patternRes <- patternTest(sce)
  earlyDETest(sce, knots = c(1, 2))
}

plotGeneCount(sds, counts, gene = sigGeneStart, models = sce)
# plotSmoothers(sce, counts, gene = sigGeneStart)
topgenes <- dplyr::filter(assoRes, pvalue < 0.05) %>% dplyr::slice_max(order_by = waldStat, n = 300) %>% rownames()
time <- slingPseudotime(sds) %>% tibble::as_tibble(rownames = "cells") %>% tibble::deframe()
heatmap_trajectory(t(counts[topgenes,]), time)


## THIS IS THE MASTER SCRIPT FOR 10X SCRNASEQ PREPROCESSING FROM COUNT MATRIX
## TO A NORMALIZED (EVEN COVARIATES-REGRESSED) DATASET SAVED AS A RICH SEURAT
## OBJECT.
## ---
## This recommended workflow comes into 3 to 4 steps :
## 1) A completely unfiltered first step, just to display the first space of ALL
##    identified cells, despite putatively many bias source present in the data.
## 2) A filtered step, except for identified doublets : just to visualize them.
## 3) A filtered + doublets removed step, but with "basic" normalization (ie :
##    no covariate regressed)
## 4) Depending on the necessity, same as 3) but with one or multiple covariates 
##    regressed.
## ---
## The ACTUAL, POOR script structure contains :
## 1) A series of blocks of variables, which are run/project -specific (a
##    detailed example is provided in the first blocl, "P30-MAMA1")
## 2) A block of common variables (variables that are very rarely to get
##    modified, paths to resources, etc...)
## 3) The functions, where all the magic happens.
## 4) A series of code blocks, corresponding to examples of each of the 4
##    described steps.
## ---
## 
## ---
## WARNINGS :
## 1) PLEASE ONLY CONSIDER FUNCTIONS AS REAL CODE. CODE BLOCKS ARE LOSELY
##    MAINTAINED AND SHOULD BE BASIC HELPERS FOR THE ORDER OF FUNCTION CALLS.
## 2) The "clustering.eval.mt" function is a multithreaded version of the
##    "clustering.eval" function. Everything was done to limit its memory
##    hunger (a minimal Seurat object is temporarily created, only necessary
##    objects are exported), but running it on 4 threads can take up to 3x the
##    memory of the monothreaded version !
## 3) The use of scBFA/BPCA implies a variation of the pipeline and its steps.
##    Please read the corresponding function header.

## INSTALL :
## R3.6.x & rstudio installed via conda -c r
## packages installed from BioConductor via BiocManager::install :
## Seurat, scran, DropletUtils, scds, scDblFinder, scBFA, (celda)

## GUILTY AUTHOR : Bastien JOB (bastien.job@gustaveroussy.fr)




## FUNCTIONS
############

## multithreading cluster
create.parallel.instance <- function(nthreads = 1) {
  doParallel::registerDoParallel(nthreads)
  cl <- BiocParallel::DoparParam()
  return(cl)
}

## Loading data into a Seurat object
load.sc.data <- function(data.path = NULL, sample.name = NULL, assay = 'RNA', droplets.limit = 1E+05, emptydrops.fdr = 1E+03, return.matrix = FALSE, filter.replicates = TRUE, BPPARAM = BiocParallel::SerialParam(), my.seed = 1337) {
  if (file.exists(data.path) && !is.null(sample.name)) {
    message("Loading data ...")
    source.format <- ""
    if(file.exists(paste0(data.path, "counts/raw_feature_bc_matrix/matrix.mtx.gz"))) { ### Cell Ranger
      source.format <- "CellRanger"
      scmat <- Seurat::Read10X(paste0(data.path,"counts/raw_feature_bc_matrix"))
      if ('Gene Expression' %in% names (scmat)) {
        message("Keeping only gene expression")
        scmat <- scmat[['Gene Expression']]
      }
    } else if(file.exists(paste0(data.path, "/", sample.name, ".mtx"))) { ### BUStools
      source.format <- "BUStools"
      scmat <- BUSpaRse::read_count_output(dir = data.path, name = sample.name, tcc = FALSE)
    } else if (file.exists(paste0(data.path, "/quants_mat.gz"))) { ### Alevin
      source.format <- "Alevin"
      scmat <- Seurat::ReadAlevin(data.path)
    } else if (file.exists(paste0(data.path, "/", sample.name, "_counts.tsv.gz"))) { ### UMI-tools
      source.format <- "UMIt-ools"
      scmat <- read.table(file = paste0(data.path, "/", sample.name, "_counts.tsv.gz"), header = TRUE, sep = "\t", quote = "", check.names = FALSE, row.names = 1)
    } else {
      stop(paste0("No data found in [", data.path, "] (wrong path ?)"))
    }
    message(paste0("Found ", source.format, " data"))
    
    scmat <- scmat[,order(colnames(scmat))]
    
    message('Droplets matrix dimensions :')
    droplets.nb <- ncol(scmat)
    print(dim(scmat))
    
    ## Filtering duplicated cell barcodes (keeping the most populated entry)
    if (filter.replicates) {
      dup.bc <- unique(colnames(scmat)[duplicated(colnames(scmat))])
      
      if(length(dup.bc) > 0) {
        message(paste0("Found ", length(dup.bc), ' (', sprintf("%.2f", length(dup.bc) / ncol(scmat) * 100), "%) replicated cell barcodes ! Filtering ..."))
        
        require(dplyr)
        dedup.tbl <- tibble(barcode = colnames(scmat), ori.index = 1:ncol(scmat), count = DelayedArray::colSums(scmat)) %>% group_by(barcode) %>% dplyr::slice(which.max(count))
        scmat <- scmat[, dedup.tbl$ori.index]
        rm(dedup.tbl)
        
        message('Droplets matrix dimensions (deduplicated) :')
        droplets.nb <- ncol(scmat)
        print(dim(scmat))
        
      } else  message('No replicated barcode found.')
    }
    
    message('Total UMIs :')
    umi.total.nb <- sum(scmat)
    print(umi.total.nb)
    if (ncol(scmat) > droplets.limit & !is.null(emptydrops.fdr)) {
      ## Removing empty droplets
      message("Removing empty droplets with emptyDrops")
      bc_rank <- DropletUtils::barcodeRanks(scmat)
      set.seed(my.seed)
      bc_rank2 <- DropletUtils::emptyDrops(scmat, BPPARAM = BPPARAM)
      scmat <- scmat[, which(bc_rank2$FDR < emptydrops.fdr)]
      message('Droplets matrix dimensions (filtered) :')
      print(dim(scmat))
      message('Total UMIs (filtered) :')
      umi.kept.nb <- sum(scmat)
      print(umi.kept.nb)
      message('Fraction of UMIs in cells :')
      print(umi.kept.nb / umi.total.nb)
      
      ## Cleaning
      rm(bc_rank, bc_rank2)
    } else {
      umi.kept.nb <- umi.total.nb
    }
    
    if (return.matrix) return(scmat)
    
    sobj <- Seurat::CreateSeuratObject(counts = scmat, project = sample.name, assay = assay)
    rm(scmat)
    # sobj$log_nCount_RNA <- log(sobj$nCount_RNA)
    sobj[[paste0('log_nCount_', assay)]] <- log(sobj[[paste0('nCount_', assay)]])
    sobj@misc$droplets <- droplets.nb
    sobj@misc$umi.raw <- umi.total.nb
    sobj@misc$umi.filtered <- umi.kept.nb
    sobj@misc$cells.ori <- ncol(sobj)
    sobj@misc$params <- list(emptydrops.fdr = emptydrops.fdr, seed = my.seed)
    sobj@misc$Rsession <- sessionInfo()
    return(sobj)
  } else stop('Data source does not exist, or no sample name specified !')
}

## Basic QC metrics
QC.metrics <- function(sobj = NULL, mt.genes.file = NULL, crb.genes.file = NULL, str.genes.file = NULL, pcmito.range = c(0, .1), pcribo.range = c(0, 1), min.features = 200, min.counts = 1000, nbin = 10, BPPARAM = BiocParallel::SerialParam()) {
  if(!is.null(sobj)) {
    if(!is.null(mt.genes.file)) if (!file.exists(mt.genes.file)) stop('mt.genes.file not found !')
    if(!is.null(crb.genes.file)) if (!file.exists(crb.genes.file)) stop('crb.genes.file not found !')
    if(!is.null(str.genes.file)) if (!file.exists(str.genes.file)) stop('str.genes.file not found !')
    
    my.seed <- sobj@misc$params$seed
    
    pcQC <- scater::perCellQCMetrics(Seurat::as.SingleCellExperiment(sobj), BPPARAM = BPPARAM)
    sobj@meta.data <- cbind(sobj@meta.data, as.data.frame(pcQC[,grep("percent", colnames(pcQC))]))
    
    ### MITO
    if (!is.null(mt.genes.file)) {
      ## Manual percent
      mito.symbols <- readRDS(mt.genes.file)
      sobj@misc$mito.symbols <- mito.symbols
      sobj@misc$params$pcmito.range = pcmito.range
      inmito <- rownames(sobj@assays$RNA@counts) %in% mito.symbols
      sobj$percent.mt <- as.vector(Matrix::colSums(sobj@assays$RNA@counts[inmito,]) / sobj$nCount_RNA)
      message('% Mitochondrial expression, ')
      print(summary(sobj$percent.mt))
      pcmito_leftbound <- sobj$percent.mt >= pcmito.range[1]
      pcmito_rightbound <- sobj$percent.mt <= pcmito.range[2]
      pcmito_num <- as.numeric(pcmito_leftbound) + (as.numeric(pcmito_rightbound)+2)
      sobj$pcmito_inrange <- pcmito_leftbound & pcmito_rightbound
      message(paste0(pcmito.range[1], ' <= % mito <= ', pcmito.range[2], ' :'))
      pcmito_factor <- factor(as.numeric(pcmito_leftbound) + abs(as.numeric(pcmito_rightbound)-1))
      levels(pcmito_factor)[levels(pcmito_factor) == 0] <- "out.left"
      levels(pcmito_factor)[levels(pcmito_factor) == 1] <- "in"
      levels(pcmito_factor)[levels(pcmito_factor) == 2] <- "out.right"
      print(table(pcmito_factor))
      ## Seurat AddModuleScore
      sobj@meta.data$MTscore <- Seurat::AddModuleScore(object = sobj, features = list(mito.symbols), nbin = nbin, seed = my.seed, assay = 'RNA', name = 'SCORE')@meta.data$SCORE1
    }
    ### RIBO
    if (!is.null(crb.genes.file)) {
      ## Manual percent
      ribo.symbols <- readRDS(crb.genes.file)
      sobj@misc$ribo.symbols <- ribo.symbols
      sobj@misc$params$pcribo.range = pcribo.range
      inribo <- rownames(sobj@assays$RNA@counts) %in% ribo.symbols
      sobj$percent.rb <- as.vector(Matrix::colSums(sobj@assays$RNA@counts[inribo,]) / sobj$nCount_RNA)
      message('% Ribosomal expression :')
      print(summary(sobj$percent.rb))
      pcribo_leftbound <- sobj$percent.rb >= pcribo.range[1]
      pcribo_rightbound <- sobj$percent.mt <= pcribo.range[2]
      sobj$pcribo_inrange <- pcribo_leftbound & pcribo_rightbound
      message(paste0(pcribo.range[1], ' <= % ribo <= ', pcribo.range[2], ' :'))
      pcribo_factor <- factor(as.numeric(pcribo_leftbound) + abs(as.numeric(pcribo_rightbound)-1))
      levels(pcribo_factor)[levels(pcribo_factor) == 0] <- "out.left"
      levels(pcribo_factor)[levels(pcribo_factor) == 1] <- "in"
      levels(pcribo_factor)[levels(pcribo_factor) == 2] <- "out.right"
      print(table(pcribo_factor))
      ## Seurat AddModuleScore
      sobj@meta.data$RBscore <- Seurat::AddModuleScore(object = sobj, features = list(ribo.symbols), nbin = nbin, seed = my.seed, assay = 'RNA', name = 'SCORE')@meta.data$SCORE1
    }
    ### STRESS
    if (!is.null(str.genes.file)) {
      ## Manual percent
      stress.symbols <- readRDS(str.genes.file)
      sobj@misc$stress.symbols <- stress.symbols
      instress <- rownames(sobj@assays$RNA@counts) %in% stress.symbols
      sobj$percent.st <- as.vector(Matrix::colSums(sobj@assays$RNA@counts[instress,]) / sobj$nCount_RNA)
      message('% Stress response expression, ')
      print(summary(sobj$percent.st))
      ## Seurat AddModuleScore
      sobj@meta.data$STscore <- Seurat::AddModuleScore(object = sobj, features = list(stress.symbols), nbin = nbin, seed = my.seed, assay = 'RNA', name = 'SCORE')@meta.data$SCORE1
    }
    
    sobj$min_features <- sobj$nFeature_RNA >= min.features
    message(paste0('% number of features >= ', min.features, ' :'))
    print(table(sobj$min_features))
    
    sobj$min_counts <- sobj$nCount_RNA >= min.counts
    message(paste0('% total counts >= ', min.counts, ' :'))
    print(table(sobj$min_counts))

    sobj@misc$params$min.features = min.features
    sobj@misc$params$min.counts = min.counts
    sobj@misc$Rsession <- sessionInfo()
  }
  return(sobj)
}

## Plot QC histograms (after using QC.metrics()  or cells.QC.filters())
QC.hist <- function(sobj = NULL, out.dir = NULL) {
  if ("Seurat" %in% is(sobj) && !is.null(out.dir)) {
    require(patchwork)
    sample.name <- Seurat::Project(sobj)
    # dir.create(out.dir, recursive = TRUE, showWarnings = TRUE)
    # if ("RNA" %in% names(sobj@assays)) {
      histFEAT <- ggplot2::qplot(sobj$nFeature_RNA, geom = "histogram", bins = 101, fill = I("white"), col = I("black"), main=paste0("nFeature_RNA (>= ", sobj@misc$params$min.features, " : ", length(which(sobj$min_features)), " cells)"), xlab = "nFeature_RNA") + ggplot2::geom_vline(xintercept=sobj@misc$params$min.features, col = "red", linetype = "dashed", size = 2) + Seurat::DarkTheme()
      histNC <- ggplot2::qplot(sobj$nCount_RNA, geom = "histogram", bins = 101, fill = I("white"), col = I("black"), main=paste0("nCount_RNA (>= ", sobj@misc$params$min.counts, " : ", length(which(sobj$min_counts)), " cells)"), xlab = "nCount_RNA") + ggplot2::geom_vline(xintercept=sobj@misc$params$min.counts, col = "red", linetype = "dashed", size = 2) + Seurat::DarkTheme()
      histMT <- ggplot2::qplot(sobj$percent.mt, geom = "histogram", bins = 101, fill = I("white"), col = I("black"), main=paste0("%mito (in [", sobj@misc$params$pcmito.range[1], ';', sobj@misc$params$pcmito.range[2], "] : ", length(which(sobj$pcmito_inrange)), " cells)"), xlab = "%mito") + ggplot2::geom_vline(xintercept=sobj@misc$params$pcmito.range, col = "red", linetype = "dashed", size = 2) + Seurat::DarkTheme()
      histRB <- ggplot2::qplot(sobj$percent.rb, geom = "histogram", bins = 101, fill = I("white"), col = I("black"), main=paste0("%ribo (in [", sobj@misc$params$pcribo.range[1], ';', sobj@misc$params$pcribo.range[2], "] : ", length(which(sobj$pcribo_inrange)), " cells)"), xlab = "%ribo") + ggplot2::geom_vline(xintercept=sobj@misc$params$pcribo.range, col = "red", linetype = "dashed", size = 2) + Seurat::DarkTheme()
    # }
    # if ("SCT" %in% names(sobj@assays)) {
    #   histFEAT2 <- ggplot2::qplot(sobj$nFeature_SCT, geom = "histogram", bins = 101, fill = I("white"), col = I("black"), main="nFeature_SCT", xlab = "nFeature (norm)") + Seurat::DarkTheme()
    #   histNC2 <- ggplot2::qplot(sobj$nCount_SCT, geom = "histogram", bins = 101, fill = I("white"), col = I("black"), main="nCount_SCT", xlab = "nCount (norm)") + Seurat::DarkTheme()
    #   histMT2 <- ggplot2::qplot(sobj$percent.mt2, geom = "histogram", bins = 101, fill = I("white"), col = I("black"), main="%mito_SCT", xlab = "%mito (norm)") + ggplot2::geom_vline(xintercept=sobj@misc$params$pcmito.range, col = "red", linetype = "dashed", size = 2) + Seurat::DarkTheme()
    #   histRB2 <- ggplot2::qplot(sobj$percent.rb2, geom = "histogram", bins = 101, fill = I("white"), col = I("black"), main="%ribo_SCT", xlab = "%ribo (norm)") + ggplot2::geom_vline(xintercept=sobj@misc$params$pcribo.range, col = "red", linetype = "dashed", size = 2) + Seurat::DarkTheme()
    #   png(paste0(out.dir, '/', sample.name, '_QChist.png'), width = 1600, height = 900)
    #   print(Seurat::CombinePlots(plots = list(histFEAT, histNC, histMT, histRB, histFEAT2, histNC2, histMT2, histRB2), ncol = 4))
    #   dev.off()
    # } else {
      png(paste0(out.dir, '/', sample.name, '_QChist.png'), width = 1600, height = 500)
      # print(Seurat::CombinePlots(plots = list(histFEAT, histNC, histMT, histRB), ncol = 4))
      print(histFEAT + histNC + histMT + histRB + plot_layout(ncol = 4))
      dev.off()
    # }
  }
}

## Filter cells based on QC metrics
cells.QC.filter <- function(sobj = NULL, min.features = 500, min.counts = 2000, pcmito.range = c(0, .1), pcribo.range = c(0, 1)) {
  if (!is.null(sobj) & all(c('nFeature_RNA', 'nCount_RNA', 'percent.mt', 'percent.rb') %in% colnames(sobj@meta.data))) {
    message('Cell expression matrix dimensions (unfiltered) :')
    print(dim(sobj))
    cellskeep <- sobj$nFeature_RNA >= min.features & sobj$nCount_RNA >= min.counts & sobj$percent.mt >= pcmito.range[1] & sobj$percent.mt <= pcmito.range[2] & sobj$percent.rb >= pcribo.range[1] & sobj$percent.rb <= pcribo.range[2]
    sobj <- sobj[, cellskeep]
    message('Cell expression matrix dimensions (filtered) :')
    print(dim(sobj))
    
    sobj@misc$params$pcmito.range = pcmito.range
    sobj@misc$params$pcribo.range = pcribo.range
    sobj@misc$params$min.features = min.features
    sobj@misc$params$min.counts = min.counts
    sobj@misc$Rsession <- sessionInfo()
    
  }
  return(sobj)
}

## Predict cell cycle phase
cell.cycle.predict <- function(sobj = NULL, assay = 'RNA', cc.pairs.file = NULL, cc.genes.file = NULL, BPPARAM = NULL, nbin = 24) {
  if(!is.null(sobj)) {
    if (!file.exists(cc.pairs.file)) stop('cc.pairs.file not found !')
    if (!file.exists(cc.genes.file)) stop('cc.genes.file not found !')
    
    if(!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist !'))
    
    my.seed <- sobj@misc$params$seed
    
    cc_pairs <- readRDS(cc.pairs.file)
    cc_genes <- readRDS(cc.genes.file)
    
    sobj <- Seurat::CellCycleScoring(object = sobj, s.features = cc_genes$s.genes, g2m.features = cc_genes$g2m.genes, assay = assay, nbin = nbin, seed = my.seed)
    sobj$Seurat.S.Score <- sobj$S.Score
    sobj$Seurat.G2M.Score <- sobj$G2M.Score
    sobj$Seurat.SmG2M.Score <- sobj$S.Score - sobj$G2M.Score
    sobj$Seurat.Phase <- sobj$Phase
    
    set.seed(my.seed)
    cycres <- scran::cyclone(Seurat::as.SingleCellExperiment(sobj, assay = assay), pairs=cc_pairs, BPPARAM = BPPARAM, verbose = TRUE)
    sobj$Phase <- cycres$phases
    sobj$nG1.Score <- cycres$normalized.scores$G1
    sobj$nS.Score <- cycres$normalized.scores$S
    sobj$nG2M.Score <- cycres$normalized.scores$G2M
    sobj$nSmG2M.Score <- cycres$normalized.scores$S - cycres$normalized.scores$G2M
    sobj$G1.Score <- cycres$scores$G1
    sobj$S.Score <- cycres$scores$S
    sobj$G2M.Score <- cycres$scores$G2M
    sobj$SmG2M.Score <- cycres$scores$S - cycres$scores$G2M
    message("Cell cycle phases : ")
    print(table(sobj$Phase))
    
    # sobj@misc$params$cc_pairs = cc_pairs;
    sobj@misc$Rsession <- sessionInfo()
  }
  return(sobj)
}

## Features filtering on fixed criteria
features.filter <- function(sobj = NULL, min.cells = 5)  {
  if(!is.null(sobj)) {
    if(all(Seurat::Assays(sobj) == 'RNA'))  {
      message('Cell expression matrix dimensions (unfiltered) :')
      print(dim(sobj))
      featureskeep <- Matrix::rowSums(sobj@assays$RNA@counts) >= min.cells
      sobj <- sobj[featureskeep,]
      message('Cell expression matrix dimensions (filtered) :')
      print(dim(sobj))
      
      sobj@misc$Rsession <- sessionInfo()
    } else {
      message("WARNING : Seurat object contained other assays than 'RNA', no filtering performed.")
    }
  }
  return(sobj)
}

## Finding cell doublets
find.doublets <- function(sobj = NULL, assay = 'RNA', min.clust.size = 100) { ## Here sobj is a SCE object
  if (!is.null(sobj)) {
    if(!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist !'))
    
    my.seed <- sobj@misc$params$seed
    
    message('Cell expression matrix dimensions :')
    print(dim(sobj))
    actual.min.clust.size <- min(min.clust.size, round(ncol(sobj) / 50))
    
    ## scran
    sobjSCE <- Seurat::as.SingleCellExperiment(sobj, assay = assay)
    set.seed(my.seed)
    qclust <- scran::quickCluster(sobjSCE, min.size = min.clust.size, use.ranks = TRUE, BPPARAM = cl)
    table(qclust)
    sobjSCE <- scran::computeSumFactors(sobjSCE, cluster = qclust, BPPARAM = cl)
    sobj$sizeFactors <- sobjSCE$sizeFactors <- DESeq2::sizeFactors(sobjSCE)
    set.seed(my.seed)
    sobj$scran.doubletscore <- scran::doubletCells(sobjSCE, size.factors.norm = sobj$sizeFactors, force.match = TRUE, BSPARAM = BiocSingular::IrlbaParam())
    rm(sobjSCE)
    sobj$log_scran.doubletscore <- log(sobj$scran.doubletscore+1)
    
    ## scDblFinder
    set.seed(my.seed)
    sobj$scDblFinder.class <- Seurat::as.Seurat(scDblFinder::scDblFinder(Seurat::as.SingleCellExperiment(sobj, assay = assay), minClusSize = actual.min.clust.size, dbr = round(ncol(sobj)^2 / 1E+05)))$scDblFinder.class ## Do not try using BPPARAM, it only works when computing on multiple samples
    sobj$scDblFinder.class <- unname(sobj$scDblFinder.class == "doublet")
    message('scDblFinder doublets :')
    print(table(sobj$scDblFinder.class))
    
    ## scds
    set.seed(my.seed)
    sobj$hybrid_score <- scds::cxds_bcds_hybrid(Seurat::as.SingleCellExperiment(sobj, assay = assay))$hybrid_score
    sobj$hybrid_score.class <- unname(sobj$hybrid_score > 1)
    message('scds-hybrid doublets :')
    print(table(sobj$hybrid_score.class))
    sobj$doublets_consensus.class <- sobj$scDblFinder.class | sobj$hybrid_score.class
    message('Consensus doublets :')
    print(table(sobj$doublets_consensus.class))
    
    sobj@misc$doublets <- list(scDblFinder = length(which(sobj$scDblFinder.class)), scds_hybrid = length(which(sobj$hybrid_score.class)), union = length(which(sobj$doublets_consensus.class)))
    sobj@misc$Rsession <- sessionInfo()
  }
  return(sobj)
}

## Filtering doublets
filter.doublets <- function(sobj = NULL, method = "both") { ## Method can be 'both', 'scDblFinder', 'scds'
  if (!is.null(sobj) && method %in% c('both', 'scDblFinder', 'scds')) {
    message('Cell expression matrix dimensions (unfiltered) :')
    print(dim(sobj))
    if (method == 'both') sobj <- sobj[, !sobj$doublets_consensus.class]
    if (method == 'scDblFinder') sobj <- sobj[, !sobj$scDblFinder.class]
    if (method == 'scds') sobj <- sobj[, !sobj$hybrid_score.class]
    message('Cell expression matrix dimensions (filtered) :')
    print(dim(sobj))
    
    sobj@misc$Rsession <- sessionInfo()
  }
  return(sobj)
}

## Control genes (if any)
tag.ctrl.genes <- function(sobj = NULL, ctrl.genes = c("GAPDH"), ctrl.min.counts = 3, assay = "SCT") {
  if (!is.null(ctrl.genes) && !is.null(sobj)) {
    if(!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist !'))
    if (length(ctrl.genes) > 0) {
      ok.ctrl <- ctrl.genes %in% rownames(sobj@assays[[assay]]@counts)
      if (any(ok.ctrl)) for (ctrlg in ctrl.genes[ctrl.genes %in% rownames(sobj@assays[[assay]]@counts)]) suppressMessages(sobj[[paste0('ctrl_', assay, '_', ctrlg)]] <- sobj@assays[[assay]]@data[rownames(sobj@assays[[assay]]@data) == ctrlg] >= ctrl.min.counts)
      
      sobj@misc$params$ctrl.genes <- ctrl.genes
      sobj@misc$params$ctrl.min.counts <- ctrl.min.counts
      sobj@misc$Rsession <- sessionInfo()
    }
  }
  return(sobj)
}

## Normalization (Seurat SCTransform)
sc.normalization <- function(sobj = NULL, assay = 'RNA', method = "SCTransform", features.n = 3000, vtr = NULL) {  ## Here sobj is a Seurat object
  if (!is.null(sobj)) {
    
    my.seed <- sobj@misc$params$seed
    
    if(!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist !'))
    if(!is.null(vtr)) vtr <- sort(vtr)
    if (toupper(method) == toupper("SCTransform")) {
      sobj <- Seurat::SCTransform(object = sobj, assay = assay, seed.use = my.seed, variable.features.n = features.n, vars.to.regress = vtr, return.only.var.genes = TRUE)
      # sobj@assays$SCT@misc$scale.data.full <- Seurat::GetResidual(sobj, features = rownames(sobj@assays$SCT@data), clip.range = sobj@assays$SCT@misc$vst.out$arguments$res_clip_range, replace.value = TRUE)@assays$SCT@scale.data
      # inmito2 <- rownames(sobj@assays$SCT@counts) %in% sobj@misc$mito.symbols
      # inribo2 <- rownames(sobj@assays$SCT@counts) %in% sobj@misc$ribo.symbols
      # sobj$percent.mt2 <- as.vector(Matrix::colSums(sobj@assays$SCT@counts[inmito2,]) / sobj$nCount_SCT)
      # sobj$percent.rb2 <- as.vector(Matrix::colSums(sobj@assays$SCT@counts[inribo2,]) / sobj$nCount_SCT)
      sobj@misc$params$normalization[[assay]] <- list(method = method, assay.out = 'SCT', features.used = features.n)
    } else if (toupper(method) == toupper("LogNormalize")) {
      sobj <- Seurat::NormalizeData(sobj, normalization.method = 'LogNormalize', assay = assay)
      sobj <- Seurat::FindVariableFeatures(sobj, assay = assay, nfeatures = features.n)
      # sobj <- Seurat::ScaleData(sobj, assay = assay, features = rownames(sobj))
      # sobj@assays[[assay]]@misc$scale.data.full <- sobj@assays[[assay]]@scale.data
      # inmito2 <- rownames(sobj@assays$RNA@counts) %in% sobj@misc$mito.symbols
      # inribo2 <- rownames(sobj@assays$RNA@counts) %in% sobj@misc$ribo.symbols
      # sobj$percent.mt2 <- as.vector(Matrix::colSums(sobj@assays$RNA@counts[inmito2,]) / sobj$nCount_RNA)
      # sobj$percent.rb2 <- as.vector(Matrix::colSums(sobj@assays$RNA@counts[inribo2,]) / sobj$nCount_RNA)
      sobj@misc$params$normalization[[assay]] <- list(method = method, assay.out = assay, features.used = features.n)
    } else if (toupper(method) == toupper("CLR")) {
      sobj <- Seurat::NormalizeData(sobj, normalization.method = 'CLR', assay = assay)
      sobj <- Seurat::FindVariableFeatures(sobj, assay = assay, nfeatures = features.n)
      # sobj <- Seurat::ScaleData(sobj, features = rownames(sobj))
      # sobj@assays[[assay]]@misc$scale.data.full <- sobj@assays[[assay]]@scale.data
      sobj@misc$params$normalization[[assay]] <- list(method = method, assay.out = assay, features.used = features.n)
    } else stop('Unknown or unsupported normalization method !')
    sobj@misc$params$vtr <- ifelse(is.null(vtr), NA, vtr)
    sobj@misc$Rsession <- sessionInfo()
  }
  return(sobj)
}

## Dimensions reduction
dimensions.reduction <- function(sobj = NULL, reduction = 'pca', assay = 'SCT', max.dims = 100L) {
  if (!is.null(sobj)) {
    
    my.seed <- sobj@misc$params$seed
    
    if(!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist !'))
    
    if (sum(dim(sobj@assays[[assay]]@scale.data)) == 0) {
      sobj <- Seurat::ScaleData(object = sobj, assay = assay)
    }
    if (reduction == 'pca') {
      sobj <- Seurat::RunPCA(object = sobj, assay = assay, verbose = FALSE, npcs = max.dims, reduction.name = paste0(assay, '_', reduction), reduction.key = paste0(assay, reduction, '_'), seed.use = my.seed)
    } else if (reduction == 'bpca') {
      # sobj <- Seurat::GetResidual(object = sobj, features = rownames(sobj@assays[[assay]]@data), assay = assay, verbose = FALSE)
      set.seed(my.seed)
      bpca.res <- scBFA::BinaryPCA(Seurat::as.SingleCellExperiment(sobj, assay = assay), X = NULL)
      sobj@reductions[[paste0(assay, '_', reduction)]] <- Seurat::CreateDimReducObject(embeddings = bpca.res$x, loadings = bpca.res$rotation, assay = assay, stdev = bpca.res$sdev, key = paste0(assay, reduction, '_'), misc = list(center = bpca.res$center, scale = bpca.res$scale))
      rm(bpca.res)
    } else if (reduction == 'bfa') { ## WARNING : DID NOT TEST IF scBFA REALLY TAKES SCT object scale.data or data !!
      set.seed(my.seed)
      bfa.res <- scBFA::scBFA(scData = as.matrix(sobj@assays[[assay]]@counts), numFactors = max.dims, X = NULL)
      dimnames(bfa.res$ZZ) <- list(colnames(sobj@assays[[assay]]@counts), paste0(paste0(assay, '_', reduction, '_'), 1L:max.dims))
      dimnames(bfa.res$AA) <- list(rownames(sobj@assays[[assay]]@counts), paste0(paste0(assay, '_', reduction, '_'), 1L:max.dims))
      sobj@reductions[[paste0(assay, '_', reduction)]] <- Seurat::CreateDimReducObject(embeddings = bfa.res$ZZ, loadings = bfa.res$AA, assay = assay, stdev = matrixStats::colSds(bfa.res$ZZ), key = paste0(assay, reduction, '_'), misc = list())
      rm(bfa.res)
    } else if (reduction == 'ica') {
      sobj <- Seurat::RunICA(object = sobj, assay = assay, verbose = FALSE, nics = max.dims, reduction.name = paste0(assay, '_', reduction), reduction.key = paste0(assay, reduction, '_'), seed.use = my.seed)
    } else if (reduction == 'mds') {
      if (assay == 'SCT') sobj <- Seurat::GetResidual(object = sobj, features = rownames(sobj@assays[[assay]]@data), assay = assay, verbose = FALSE)
      set.seed(my.seed)
      # mds.res <- scater::calculateMDS(x = Seurat::as.SingleCellExperiment(x = sobj, assay = assay), ncomponents = max(red.dims))
      mds.res <- scater::calculateMDS(x = sobj@assays[[assay]]@scale.data, ncomponents = max.dims)
      sobj@reductions[[paste0(assay, '_', reduction)]] <- Seurat::CreateDimReducObject(embeddings = mds.res, loadings = matrix(nrow = 0, ncol = 0), assay = assay, stdev = matrixStats::colSds(mds.res), key = paste0(assay, reduction, '_'), misc = list())
    } else stop("Unknown reduction method !")
    sobj@assays[[assay]]@scale.data <- matrix(nrow = 0, ncol = 0)
  }
  return(sobj)
}

## Special function for scBFA / BPCA processing (replaces sc.normalization(), dimensions.reduction()).
## Performs LogNormalization or SCTransform (norm parameter)
## 1) LogNormalization or SCTransform is performed to reduces libsize effect
## 2) HVG detection is performed for LogNorm (included in SCtransform) : BFA/BPCA is more efficient on HVGs
## 3) reduction is performed on matrix restrained to HVGs
## 'method' should be : 'scbfa' or 'bpca'.
## 'X' is a (cell x score) matrix of variables to regress
## 'max.dims is for scBFA which is way slower, so we can cut at a defined number of dimensions
binary.processing <- function(sobj = NULL, assay = 'RNA', norm.method = 'LogNormalize', method = 'scbfa', features.n = 3000, max.dims = 30, vtr = NULL, vtr.scale = TRUE) {
  if (is.null(sobj)) stop('No Seurat object provided !')
  if(!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist !'))
  
  if (!(tolower(method) %in% c('scbfa', 'bpca'))) stop("Requested method is not supported !")
  
  my.seed <- sobj@misc$params$seed
  
  ## Checking VTR
  if (!is.null(vtr)) {
    if (!all(vtr %in% colnames(sobj@meta.data))) stop('Not all vtr names found in object meta data !')
    X <- sobj@meta.data[,colnames(sobj@meta.data) %in% vtr, drop = FALSE]
    for(x in seq_len(ncol(X))) {
      if (!is.numeric(X[,x])) X[,x] <- as.numeric(as.factor(X[,x]))
    }
    X <- as.matrix(X)
    # X <- as.matrix(sobj@meta.data[, colnames(sobj@meta.data) %in% vtr])
  } else X <- vtr
  
  if (tolower(norm.method) == 'lognormalize') {
    ## LN
    # sobj[[paste0(assay, "_RAW")]] <- Seurat::CreateAssayObject(counts = sobj@assays[[assay]]@counts)
    sobj <- Seurat::NormalizeData(sobj, normalization.method = 'LogNormalize', assay = assay)
    
    ## HVG
    sobj <- Seurat::FindVariableFeatures(sobj, assay = assay, nfeatures = features.n)
    
  } else if (tolower(norm.method) == 'sctransform') {
    sobj <- Seurat::SCTransform(sobj, assay = assay, variable.features.n = features.n, seed.use = my.seed)
    assay <- 'SCT'
  } else stop("Unsupported normalization method !")
  
  ## Scaling
  # sobj <- Seurat::ScaleData(sobj, assay = assay, features = rownames(sobj))
  # sobj@assays[[assay]]@misc$scale.data.full <- sobj@assays[[assay]]@scale.data
  
  ## extracting HVGs index
  hvg.idx <- rownames(sobj@assays[[assay]]@counts) %in% sobj@assays[[assay]]@var.features
  
  ## Scaling vtr if requested
  if (!is.null(X) & vtr.scale) X <- scale(X)
  
  ## Reduction
  my.red <- paste0(assay, '_', reduction)
  if (tolower(method) == "scbfa") {
    set.seed(my.seed)
    bfa.res <- scBFA::scBFA(scData = as.matrix(sobj@assays[[assay]]@counts[hvg.idx,]), numFactors = max.dims, X = X)
    dimnames(bfa.res$ZZ) <- list(colnames(sobj@assays[[assay]]@counts), paste0('SCBFA_', 1L:max.dims))
    dimnames(bfa.res$AA) <- list(rownames(sobj@assays[[assay]]@counts)[hvg.idx], paste0('SCBFA_', 1L:max.dims))
    sobj@reductions[[my.red]] <- Seurat::CreateDimReducObject(embeddings = bfa.res$ZZ, loadings = bfa.res$AA, assay = assay, stdev = matrixStats::colSds(bfa.res$ZZ), key = paste0(assay, reduction, '_'), misc = list())
    sobj@misc$params$reduction[[my.red]] <- list(numFactors = bfa.res$numFactors, X = if(is.null(X)) NA else colnames(X))
    # sobj@reductions$scbfa@misc$from.assay <- assay ## Already in sobj@reductions$scbfa@assay.used
    sobj@reductions[[my.red]]@misc$binary.matrix <- bfa.res$BB
    rm(bfa.res)
    
  } else if (tolower(method) == "bpca") {
    set.seed(my.seed)
    bpca.res <- scBFA::BinaryPCA(scData = as.matrix(sobj@assays[[assay]]@counts[hvg.idx,]), X = X)
    colnames(bpca.res$x) <- paste0('BPCA_', 1L:features.n)
    colnames(bpca.res$rotation) <- paste0('BPCA_', 1L:features.n)
    sobj@reductions[[my.red]] <- Seurat::CreateDimReducObject(embeddings = bpca.res$x, loadings = bpca.res$rotation, assay = assay, stdev = bpca.res$sdev, key = paste0(assay, reduction, '_'), misc = list(center = bpca.res$center, scale = bpca.res$scale))
    sobj@misc$params$reduction[[my.red]] <- list(numFactors = bpca.res$numFactors, X = if(is.null(X)) NA else colnames(X))
    # sobj@reductions$bpca@misc$from.assay <- assay ## Already in sobj@reductions$scbfa@assay.used
    rm(bpca.res)
    
  } else stop(paste0('Method "', method, '" is not supported !'))
  
  rm(hvg.idx)
  
  ## Misc
  # sobj@assays[[assay]]@misc$scale.data.full <- sobj@assays[[assay]]@scale.data
  sobj@misc$params$normalization[[assay]] <- list(method = tolower(norm.method), assay.out = assay, features.used = features.n)
  sobj@reductions[[my.red]]@misc$vtr <- if (!is.null(vtr)) vtr else NA

  return(sobj)
}

## DecontX processing
### This requires a seurat object with clusters (see celda::decontX)
### This function returns a new seurat object with corrected raw count matrix
decontx.process <- function(sobj, assay = 'RNA', idents = NULL, ...) {
  if (is.null(sobj)) stop("A Seurat object is required !")
  if (is.null(idents)) stop("An ident name is required ! Use get.idents() to get available idents.")
  my.counts = as.matrix(sobj@assays[[assay]]@counts)
  mode(my.counts) = "integer"
  my.z = as.factor(as.character(as.numeric(sobj@meta.data[[idents]])))
  dx.res <- celda::decontX(counts = my.counts, z = my.z, seed = sobj@misc$params$seed, ...)
  sobj@assays[[assay]]@counts <- as(dx.res$resList$estNativeCounts, "dgCMatrix")
  rm(my.counts, dx.res, my.z)
  return(sobj)
}

## Reduction dims correlation with bias sources
dimensions.eval <- function(sobj = NULL, assay = 'RNA', reduction = 'pca', cor.method = 'spearman', meta.names = c('sizeFactors', 'nFeature_RNA', 'percent.mt', 'MTscore', 'percent.rb', 'RBscore', 'percent.st', 'STscore', "S.Score", "G1.Score", "G2M.Score", "SmG2M.Score"), eval.markers = c('GAPDH'), max.dims = 100L, out.dir = NULL) {
  if (is.null(out.dir)) stop('No output dir provided !')
  if (!dir.exists(out.dir)) stop('Output directory does not exist !')
  if (is.null(sobj)) stop('No Seurat object provided !')
  my.red <- paste0(assay, '_', reduction)
  if (!(my.red %in% names(sobj@reductions))) stop(paste0('Reduction "', my.red, '" not present in the provided Seurat object !'))
  if(!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist !'))
  
  sample.name <- Seurat::Project(sobj)
  
  ndims <- min(max.dims, ncol(sobj@reductions[[my.red]]@cell.embeddings))
  if (!dir.exists(out.dir)) dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
  if (assay == 'SCT') {
    ## Regenerate full scaled matrix
    scale.data.full <- Seurat::GetResidual(sobj, features = rownames(sobj@assays$SCT@data), clip.range = sobj@assays$SCT@misc$vst.out$arguments$res_clip_range, replace.value = TRUE)@assays$SCT@scale.data
    eval.markers <- eval.markers[eval.markers %in% rownames(scale.data.full)]
  }
  if (assay == 'RNA') eval.markers <- eval.markers[eval.markers %in% rownames(sobj@assays[[assay]]@data)]
  # print(eval.markers)
  png(paste0(out.dir, '/', sample.name, '_', reduction, '_dims.bias.cor.png'), width = 1400, height = 800)
  suppressWarnings(plot(0, 0, log = 'x', xlim = c(1L, ndims), ylim = c(0,1), type = "n", xlab = paste0(reduction, " dimension"), ylab = paste0(cor.method, " correlation"), xaxs = "i", yaxs = "i"))
  meta.names <- meta.names[meta.names %in% colnames(sobj@meta.data)]
  meta.cols <- if (length(meta.names) > 12) scales::hue_pal()(length(meta.names)) else RColorBrewer::brewer.pal(length(meta.names), name = "Paired")
  
  for(myf in seq_along(meta.names)) {
    corvec <- vapply(seq_len(ndims), function(k) { (abs(cor(sobj@reductions[[my.red]]@cell.embeddings[,k], sobj[[meta.names[myf]]], method = cor.method))) }, .1)
    lines(corvec, type = "l", ylim = c(0,1), col = meta.cols[myf], lwd = 4)
  }
  if (!is.null(eval.markers)) {
    for(mym in seq_along(eval.markers)) {
      if (assay == 'SCT') corvec <- vapply(seq_len(ndims), function(k) { (abs(cor(sobj@reductions[[my.red]]@cell.embeddings[,k], scale.data.full[rownames(sobj@assays[[assay]]@data) == eval.markers[mym],], method = cor.method))) }, .1)
      if (assay == 'RNA') corvec <- vapply(seq_len(ndims), function(k) { (abs(cor(sobj@reductions[[my.red]]@cell.embeddings[,k], sobj@assays[[assay]]@data[rownames(sobj@assays[[assay]]@data) == eval.markers[mym],], method = cor.method))) }, .1)
      
      # corvec <- vapply(seq_len(ndims), function(k) { (abs(cor(sobj@reductions[[reduction]]@cell.embeddings[,k], sobj@assays[[assay]]@data[rownames(sobj@assays[[assay]]@data) == eval.markers[mym],], method = cor.method))) }, .1)
      lines(corvec, type = "l", ylim = c(0,1), col = 1, lwd = 5, lty = mym+1L)
      # }
    }
  }
  
  # abline(h = 0, lty = 2)
  legend(x = ndims, y = 1, legend = c(meta.names, eval.markers), col = c(meta.cols, rep(1, length(eval.markers))), lty = c(rep(1, length(meta.names)), seq_along(eval.markers)), lwd = 5, xjust = 1, yjust = 1)
  dev.off()
}

## Evaluating Louvain clusters (pca dims / resolution) (multithreading version)
clustering.eval.mt <- function(sobj = NULL, reduction = 'pca', assay = 'SCT', dimsvec = seq.int(10L,120L,10L), resvec = seq(.1,1,.1), BPPARAM = BiocParallel::SerialParam(), out.dir = NULL, solo.pt.size = 2) {
  if (is.null(out.dir)) stop('No output dir provided !')
  if (!dir.exists(out.dir)) stop('Output directory does not exist !')
  if (is.null(sobj)) stop('No Seurat object provided !')
  if (!"seed" %in% names(sobj@misc$params)) stop("No seed found in Seurat object !")
  
  my.seed <- sobj@misc$params$seed
  
  if (!(paste0(assay, '_', reduction) %in% names(sobj@reductions))) stop(paste0('Reduction "', paste0(assay, '_', reduction), '" not present in the provided Seurat object !'))
  if(!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist !'))
  if (max(dimsvec) > ncol(sobj@reductions[[paste0(assay, '_', reduction)]]@cell.embeddings)) stop(paste0('Max dimsvec requested is ', max(dimsvec), ' whereas max dimension in "', paste0(assay, '_', reduction), '" reduction is ', ncol(sobj@reductions[[paste0(assay, '_', reduction)]]@cell.embeddings)))
  
  sample.name <- Seurat::Project(sobj)
  clustree.dir <- paste0(out.dir, "/clustree_", reduction, '/')
  umaps.clustree.dir <- paste0(clustree.dir, "/uMAPs/")
  pca.clustree.dir <- paste0(clustree.dir, "/dimensions/")
  res.clustree.dir <- paste0(clustree.dir, "/louvain_resolution/")
  dir.create(umaps.clustree.dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(pca.clustree.dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(res.clustree.dir, recursive = TRUE, showWarnings = FALSE)
  require(ggplot2)
  require(clustree)
  
  plot.pix <- 600
  
  `%dopar%` <- foreach::"%dopar%"
  `%mydo%` <- if (BiocParallel::bpworkers(BPPARAM) > 1) foreach::"%dopar%" else foreach::"%do%"
  `%do%` <- foreach::"%do%"
  
  ## Builiding minimal Seurat object (for memory sake)
  miniobj <- sobj
  ### Removing other assays
  Seurat::DefaultAssay(miniobj) <- assay
  other.assays <- names(miniobj@assays)[names(miniobj@assays) != assay]
  if (length(other.assays) > 0) miniobj@assays[other.assays] <- NULL
  ### Removing scale.data slot
  miniobj@assays[[assay]]@scale.data <- matrix(nrow = 0, ncol = 0)
  ### Removing misc
  miniobj@assays[[assay]]@misc <- list()
  ### Removing other reductions
  other.reductions <- names(miniobj@reductions)[names(miniobj@reductions) != paste0(assay, '_', reduction)]
  if (length(other.reductions) > 0) miniobj@reductions[other.reductions] <- NULL
  ### Slimming other slots
  miniobj@commands <- miniobj@misc <- list()
  require(Matrix)
  miniobj@assays[[assay]]@counts <- new("dgCMatrix")
  
  ## Purging putative residues from a former run
  miniobj@meta.data[, grep(colnames(miniobj@meta.data), pattern = "seurat_clusters_LE_")] <- NULL
  rm(sobj)
  gc()
  
  ## Reordering dimsvec if needed (as higher dims are longer to compute)
  if (which.max(dimsvec) != 1) dimsvec <- sort(dimsvec, decreasing = TRUE)
  
  my.red <- paste0(assay, '_', reduction)
  
  # for (my.dims in dimsvec) {
  # resclust.all <- foreach::foreach(my.dims = dimsvec, .combine = "cbind") %do% {
  resclust.all <- foreach::foreach(my.dims = dimsvec, .combine = "cbind", .inorder = FALSE, .errorhandling = "stop", .packages = c("Seurat", "ggplot2")) %mydo% {
    
    miniobj@graphs <- list()
    message(paste0("Dimensions 1 to ", my.dims))
    suppressMessages(miniobj <- Seurat::FindNeighbors(object = miniobj, assay = assay, dims = 1L:my.dims, reduction = my.red))
    
    # print(resvec)
    resloop = list()
    # resloop <- foreach::foreach(my.res = resvec, .inorder = TRUE, .errorhandling = "stop", .noexport = objects(), export = c('miniobj', 'assay', 'my.seed', 'reduction', 'my.dims', 'resvec', 'umaps.clustree.dir', 'sample.name', 'solo.pt.size'), .packages = c("Seurat", "ggplot2")) %do% {
    resloop <- foreach::foreach(my.res = resvec, .inorder = FALSE, .errorhandling = "stop", .noexport = objects()) %do% {
      # for (my.res in resvec) {
      
      # print(resvec)
      message(paste0("Testing resolution ", my.res, " ..."))
      
      miniobj <- Seurat::FindClusters(object = miniobj, assay = assay, random.seed = my.seed, resolution = my.res, graph.name = paste0(assay, '_snn'))

      miniobj <- Seurat::RunUMAP(object = miniobj, assay = assay, dims = 1L:my.dims, reduction = my.red, seed.use = my.seed, reduction.name = paste(c(assay, reduction, my.dims, 'umap'), collapse = '_'))
      png(paste0(umaps.clustree.dir, '/', sample.name,'_uMAP_', reduction, my.dims, "_res", my.res*10, '.png'), width = 1100, height = 1000)
      resdim.plot <- Seurat::LabelClusters(plot = Seurat::DimPlot(object = miniobj, reduction = paste(c(assay, reduction, my.dims, 'umap'), collapse = '_'), pt.size = solo.pt.size) + ggplot2::ggtitle(paste0(toupper(reduction), " dims =  ", my.dims, " ; resolution = ", my.res)) + Seurat::DarkTheme(), id = "ident", size = solo.pt.size*3, repel = FALSE, color = "white", fontface = "bold")
      print(resdim.plot)
      dev.off()
      
      return(list(resplot = resdim.plot, clusters = miniobj$seurat_clusters))
      # resloop = c(resloop, list(resplot = resdim.plot, clusters = miniobj$seurat_clusters))
    }
    
    ## Decomposing resloop
    resplots <- unlist(resloop, recursive = FALSE)[seq.int(1, length(resloop)*2, 2)]
    resclust <- as.data.frame(unlist(resloop, recursive = FALSE)[seq.int(2, length(resloop)*2, 2)])
    resclust <- cbind(resclust, resclust)
    rm(resloop)
    colnames(resclust) <- c(paste0(paste0("seurat_clusters_LE_", reduction, my.dims, "_res", resvec)), paste0(paste0("seurat_clusters_LE_res", resvec, "_", reduction, my.dims)))
    miniobj@meta.data <- cbind(miniobj@meta.data, resclust)
    # rm(resclust)
    
    grid.xy <- grid.scalers(length(resplots))
    # png(paste0(umaps.clustree.dir, '/', sample.name, '_', assay, '_uMAPs_', reduction, my.dims, '_ALLres.png'), width = 2000, height = 1000)
    png(paste0(umaps.clustree.dir, '/', sample.name, '_', assay, '_uMAPs_', reduction, my.dims, '_ALLres.png'), width = grid.xy[1]*plot.pix, height = grid.xy[2]*plot.pix)
    print(patchwork::wrap_plots(resplots))
    dev.off()
    
    rm(resplots)
    
    if (length(resvec) > 1) {
      cres <- clustree::clustree(miniobj, prefix = paste0("seurat_clusters_LE_", reduction, my.dims, "_res"))
      png(paste0(pca.clustree.dir, '/', sample.name, '_', assay, '_', reduction, my.dims, '.png'), width = 800, height = 1000)
      print(cres + ggplot2::ggtitle(paste0(sample.name, ', ', toupper(reduction), ' = ', my.dims)))
      dev.off()
    }
    
    gc(verbose = FALSE)
    
    return(resclust)
  }
  miniobj@meta.data <- cbind(miniobj@meta.data, resclust.all)
  rm(resclust.all)
  
  for (my.res in resvec) {
    if (length(dimsvec) > 1) {
      # png(paste0(umaps.clustree.dir, '/', sample.name, '_', assay, '_uMAPS_res', my.res*10, '_ALLdims.png'), width = 2000, height = 1000)
      # print(Seurat::CombinePlots(plots = dimresplots[grep(names(dimresplots), pattern = paste0('_res', my.res))]))
      # dev.off()
      
      cres <- clustree::clustree(miniobj, prefix = paste0("seurat_clusters_LE_res", my.res, "_", reduction))
      
      png(paste0(res.clustree.dir, '/', sample.name, '_', assay, '_res', my.res*10, '.png'), width = 800, height = 1000)
      print(cres + ggplot2::ggtitle(paste0(sample.name, ", res = ", my.res)))
      dev.off()
    }
  }
  # rm(dimresplots)
  # for(x in grep(colnames(sobj@meta.data), pattern = "seurat_clusters_LE_", value = TRUE)) sobj@meta.data[[x]] <- NULL
  # for(x in grep(colnames(sobj@meta.data), pattern = "SCT_snn_res", value = TRUE)) sobj@meta.data[[x]] <- NULL
  
  # dimensions.eval(sobj, reduction = reduction, red.dims = 1L:ncol(sobj@reductions[[reduction]]), markers = eval.markers, out.dir = paste0(data.path, "/DEFAULT/"), vtr = "DEFAULT")
  rm(miniobj)
  gc(verbose = FALSE)
}

## Louvain clustering + UMAP
louvain.cluster <- function(sobj = NULL, reduction = 'pca', assay = "SCT", max.dim = 100L, resolution = .8, out.dir = NULL, solo.pt.size = 2) {
  if (is.null(out.dir)) stop('No output dir provided !')
  if (!dir.exists(out.dir)) stop('Output directory does not exist !')
  if (is.null(sobj)) stop('No Seurat object provided !')
  if (!"seed" %in% names(sobj@misc$params)) stop("No seed found in Seurat object !")
  my.red <- paste0(assay, '_', reduction)
  if (!(my.red %in% names(sobj@reductions))) stop(paste0('Reduction "', my.red, '" not present in the provided Seurat object !'))
  if(!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist !'))
  if (max.dim > ncol(sobj@reductions[[my.red]]@cell.embeddings)) stop(paste0('Max dimension requested is ', max.dim, ' whereas max dimension in "', reduction, '" reduction is ', ncol(sobj@reductions[[my.red]]@cell.embeddings)))
  
  my.seed <- sobj@misc$params$seed
    
  sample.name <- Seurat::Project(sobj)
  sobj@graphs <- list()
  set.seed(my.seed)
  sobj <- Seurat::FindNeighbors(object = sobj, assay = assay, reduction = my.red, compute.SNN = TRUE, dims = 1L:max.dim)
  # names(sobj@graphs)[names(sobj@graphs) == paste0(assay, '_nn')] <- paste0(assay, '.', reduction, '_nn')
  # names(sobj@graphs)[names(sobj@graphs) == paste0(assay, '_snn')] <- paste0(assay, '.', reduction, '_snn')
  names(sobj@graphs)[names(sobj@graphs) == paste0(assay, '_nn')] <- paste0(paste(c(assay, reduction, max.dim), collapse = '.'), '_nn')
  names(sobj@graphs)[names(sobj@graphs) == paste0(assay, '_snn')] <- paste0(paste(c(assay, reduction, max.dim), collapse = '.'))
  ori.assay <- Seurat::DefaultAssay(sobj)
  Seurat::DefaultAssay(sobj) <- assay
  # sobj <- Seurat::FindClusters(object = sobj, resolution = resolution, random.seed = my.seed, graph.name = paste0(assay, '.', reduction, '_snn'))
  sobj <- Seurat::FindClusters(object = sobj, resolution = resolution, random.seed = my.seed, graph.name = paste0(paste(c(assay, reduction, max.dim), collapse = '.')))
  Seurat::DefaultAssay(sobj) <- ori.assay
  sobj <- Seurat::RunUMAP(object = sobj, assay = assay, dims = 1L:max.dim, reduction = paste0(assay, '_', reduction), graph.name = paste0(paste(c(assay, reduction, max.dim), collapse = '.')), reduction.name = paste(c(assay, reduction, max.dim, 'umap'), collapse = "_"), reduction.key = tolower(paste0(assay, reduction, max.dim, 'umap_')), seed.use = my.seed)
  # sobj@reductions[[paste0('umap_', reduction, max.dim)]] <- sobj@reductions$umap
  # sobj <- Seurat::RunUMAP(object = sobj, assay = assay, dims = 1L:max.dim, reduction = reduction, seed.use = my.seed, n.components = 3L, reduction.name = paste(c(assay, reduction, 'umap3d', max.dim), collapse = '_'), reduction.key = paste0(assay, reduction, 'UMAP3D_'))
  sobj <- Seurat::RunUMAP(object = sobj, assay = assay, dims = 1L:max.dim, reduction = paste0(assay, '_', reduction), n.components = 3L, graph.name = paste0(assay, '.', reduction, '_snn'), reduction.name = paste(c(assay, reduction, max.dim, 'umap3d'), collapse = "_"), reduction.key = tolower(paste0(assay, reduction, max.dim, 'umap3d_')), seed.use = my.seed)
  # sobj@reductions[[paste0('umap3d_', reduction, max.dim)]] <- sobj@reductions$umap3d
  if (!is.null(out.dir)) {
    png(paste0(out.dir, '/', paste0(c(sample.name, assay, reduction), collapse = "_"), '_uMAP_dim', max.dim, "_res", resolution, '.png'), width = 1000, height = 1000)
    print(Seurat::LabelClusters(plot = Seurat::DimPlot(object = sobj, reduction = paste(c(assay, reduction, max.dim, 'umap'), collapse = "_"), pt.size = solo.pt.size) + ggplot2::ggtitle("Louvain clusters (Seurat)") + Seurat::DarkTheme(), id = "ident", size = solo.pt.size*3, repel = FALSE, color = "white", fontface = "bold"))
    dev.off()
    u3dlist <- list(
      Seurat::LabelClusters(plot = Seurat::DimPlot(object = sobj, reduction = paste(c(assay, reduction, max.dim, 'umap3d'), collapse = "_"), pt.size = solo.pt.size, dims = c(1,2)) + ggplot2::ggtitle("Louvain clusters (Seurat)") + Seurat::DarkTheme(), id = "ident", size = solo.pt.size*3, repel = FALSE, color = "white", fontface = "bold"),
      Seurat::LabelClusters(plot = Seurat::DimPlot(object = sobj, reduction = paste(c(assay, reduction, max.dim, 'umap3d'), collapse = "_"), pt.size = solo.pt.size, dims = c(1,3)) + ggplot2::ggtitle("Louvain clusters (Seurat)") + Seurat::DarkTheme(), id = "ident", size = solo.pt.size*3, repel = FALSE, color = "white", fontface = "bold"),
      Seurat::LabelClusters(plot = Seurat::DimPlot(object = sobj, reduction = paste(c(assay, reduction, max.dim, 'umap3d'), collapse = "_"), pt.size = solo.pt.size, dims = c(2,3)) + ggplot2::ggtitle("Louvain clusters (Seurat)") + Seurat::DarkTheme(), id = "ident", size = solo.pt.size*3, repel = FALSE, color = "white", fontface = "bold")
    )
    png(paste0(out.dir, '/', paste0(c(sample.name, assay, reduction), collapse = "_"), '_uMAP3D_dim', max.dim, "_res", resolution, '.png'), width = 2000, height = 2000)
    print(patchwork::wrap_plots(u3dlist) + patchwork::plot_layout(ncol = 2))
    dev.off()
    
  }
  # sobj@meta.data[[paste0("seurat_clusters_", reduction, max.dim, '_res', resolution)]] <- sobj$seurat_clusters
  sobj@misc$params$max.dim <- max.dim
  sobj@misc$params$resolution <- resolution
  sobj@misc$Rsession <- sessionInfo()
  return(sobj)
}

## Automatic annotation of cells/clusters
cells.annot <- function(sobj = NULL, assay = 'RNA', umap.name = NULL, ident = NULL, singler.setnames = NULL, clustifyr.setnames = NULL, cellassign.markers = NULL, sr.minscore = .25, cfr.minscore = .35, out.dir = NULL, solo.pt.size = 2) {
  if (is.null(out.dir)) stop('No output dir provided !')
  if (!dir.exists(out.dir)) stop('Output directory does not exist !')
  if (is.null(sobj)) stop('No Seurat object provided !')
  if (!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist !'))
  if (!"seed" %in% names(sobj@misc$params)) stop("No seed found in Seurat object !")
  if (is.null(ident)) stop("No ident name provided !")
  if (!ident %in% colnames(sobj@meta.data)) stop(paste0("Ident '", ident, "' not fond in Seurat object !"))
  if (is.null(umap.name)) stop("No umap reduction name provided !")
  if (!umap.name %in% names(sobj@reductions)) stop(paste0("UMAP reduction '", umap.name, "' not fond in Seurat object !"))
  
  ## Restoring seed and sample name from within the object
  my.seed <- sobj@misc$params$seed
  sample.name <- Seurat::Project(sobj)
  
  ## Building output structure
  cellannot.dir <- paste0(out.dir, "/cells_annotation/")
  sr.cellannot.dir <- paste0(cellannot.dir, "/singler/")
  cfr.cellannot.dir <- paste0(cellannot.dir, "/clustifyr/")
  cellassign.dir <- paste0(cellannot.dir, "/cellassign/")
  dir.create(sr.cellannot.dir, recursive = TRUE)
  dir.create(cfr.cellannot.dir, recursive = TRUE)
  dir.create(cellassign.dir, recursive = TRUE)
  
  ### SingleR
  suppressMessages(require(SingleR))
  for (setname in singler.setnames) {
    print(setname)
    singler.ref <- suppressMessages(do.call(match.fun(paste0(setname)), args = list()))
    ## per cell
    entryname <- paste0("SR_", setname, "_cells")
    singler.res <- SingleR::SingleR(test = sobj@assays[[assay]]@data, ref = singler.ref, labels = singler.ref$label.main)
    singler.res$pruned.labels[singler.res@listData$tuning.scores$first < sr.minscore] <- NA
    sobj@meta.data[[entryname]] <- as.factor(unname(singler.res$pruned.labels))
    png(paste0(sr.cellannot.dir, '/', sample.name, '_', assay, '_uMAP_', entryname, '.png'), width = 1200, height = 1200)
    if (all(is.na(sobj@meta.data[[entryname]]))) {
      suppressMessages(print(Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = solo.pt.size, cols = "grey50", group.by = entryname) + ggplot2::ggtitle(paste0("SingleR predicted cell types (", setname, ") for ", ncol(sobj), " cells")) + ggplot2::theme(legend.position = "bottom") + Seurat::DarkTheme()))
    } else {
      suppressMessages(print(Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = solo.pt.size, group.by = entryname) + ggplot2::ggtitle(paste0("SingleR predicted cell types (", setname, ") for ", ncol(sobj), " cells")) + ggplot2::theme(legend.position = "bottom") + Seurat::DarkTheme()))
    }
    dev.off()
    ## per cluster
    entryname <- paste0("SR_", setname, "_clust")
    sobj@meta.data[[entryname]] <- sobj@meta.data[[ident]]
    singler.res <- SingleR::SingleR(test = sobj@assays[[assay]]@data, ref = singler.ref, labels = singler.ref$label.main, method = "cluster", clusters = sobj[[entryname, drop = TRUE]])
    singler.res$pruned.labels[singler.res@listData$tuning.scores$first < sr.minscore] <- NA
    levels(sobj@meta.data[[entryname]]) <- singler.res$pruned.labels
    
    png(paste0(sr.cellannot.dir, '/', sample.name, '_', assay, '_uMAP_', entryname, '.png'), width = 1200, height = 1200)
    if (all(is.na(sobj@meta.data[[entryname]]))) {
      suppressMessages(print(Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = solo.pt.size, cols = "grey50", group.by = entryname) + ggplot2::ggtitle(paste0("SingleR predicted cell types (", setname, ") for ", length(unique(sobj$seurat_clusters)), " cluster(s)")) + ggplot2::theme(legend.position = "bottom") + Seurat::DarkTheme()))
    } else suppressMessages(print(Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = solo.pt.size, group.by = entryname) + ggplot2::ggtitle(paste0("SingleR predicted cell types (", setname, ") for ", nlevels(sobj@meta.data[[entryname]]), " cluster(s)")) + ggplot2::theme(legend.position = "bottom") + Seurat::DarkTheme()))
    dev.off()
  }
  
  ### CLUSTIFYR
  suppressPackageStartupMessages(require(clustifyrdata))
  for (setname in clustifyr.setnames) {
    print(setname)
    entryname <- paste0("CFR_", setname, "_cells")
    sobj@meta.data[[entryname]] <- sobj@meta.data[[ident]]
    non.var.genes <- rownames(sobj@assays[[assay]]@data)[!(rownames(sobj@assays[[assay]]@data) %in% sobj@assays[[assay]]@var.features)]
    ## per cell
    cfyr.res <- suppressMessages(clustifyr::clustify(input = sobj@assays[[assay]]@data, metadata = sobj[[entryname, drop = TRUE]], ref_mat = get0(setname), per_cell = TRUE, dr = umap.name, exclude_genes = non.var.genes))
    cfyr.res.tbl <- clustifyr::cor_to_call(cfyr.res)
    cfyr.res.tbl <- cfyr.res.tbl[order(cfyr.res.tbl$cluster),]
    sobj@meta.data[[entryname]] <- cfyr.res.tbl$type
    sobj@meta.data[[entryname]][cfyr.res.tbl$r < cfr.minscore] <- NA
    
    png(paste0(cfr.cellannot.dir, '/', sample.name, '_', assay, '_uMAP_', entryname, '.png'), width = 1200, height = 1200)
    if (all(is.na(sobj@meta.data[[entryname]]))) {
      suppressMessages(print(Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = solo.pt.size, cols = "grey50", group.by = entryname) + ggplot2::ggtitle(paste0("clustifyr predicted cell types (", setname, ") for ", ncol(sobj), " cells")) + ggplot2::theme(legend.position = "bottom") + Seurat::DarkTheme()))
    } else {
      suppressMessages(print(Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = solo.pt.size, group.by = entryname) + ggplot2::ggtitle(paste0("clustifyr predicted cell types (", setname, ") for ", ncol(sobj), " cells")) + ggplot2::theme(legend.position = "bottom") + Seurat::DarkTheme()))
    }
    dev.off()
    ## per cluster
    entryname <- paste0("CFR_", setname, "_clust")
    sobj@meta.data[[entryname]] <- sobj@meta.data[[ident]]
    cfyr.res <- suppressMessages(clustifyr::clustify(input = sobj@assays[[assay]]@data, metadata = sobj[[entryname, drop = TRUE]], ref_mat = get0(setname), per_cell = FALSE, dr = umap.name, exclude_genes = non.var.genes))
    cfyr.res.tbl <- clustifyr::cor_to_call(cfyr.res)
    cfyr.res.tbl <- cfyr.res.tbl[order(as.numeric(cfyr.res.tbl$cluster)),]
    levels(sobj@meta.data[[entryname]]) <- cfyr.res.tbl$type
    sobj@meta.data[[entryname]][cfyr.res.tbl$r < cfr.minscore] <- NA
    
    png(paste0(cfr.cellannot.dir, '/', sample.name, '_', assay, '_uMAPs_', entryname, '.png'), width = 1200, height = 1200)
    if (all(is.na(sobj@meta.data[[entryname]]))) {
      suppressMessages(print(Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = solo.pt.size, cols = "grey50", group.by = entryname) + ggplot2::ggtitle(paste0("clustifyr predicted cell types (", setname, ") for ", nlevels(sobj@meta.data[[entryname]]), " cluster(s)")) + ggplot2::theme(legend.position = "bottom") + Seurat::DarkTheme()))
    } else {
      suppressMessages(print(Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = solo.pt.size, group.by = entryname) + ggplot2::ggtitle(paste0("clustifyr predicted cell types (", setname, ") for ", nlevels(sobj@meta.data[[entryname]]), " cluster(s)")) + ggplot2::theme(legend.position = "bottom") + Seurat::DarkTheme()))
    }
    dev.off()
    
    
    
    sobj@misc$params$sr.minscore <- sr.minscore
    sobj@misc$params$cfr.minscore <- cfr.minscore
    sobj@misc$Rsession <- sessionInfo()
  }
  if (!is.null(cellassign.markers)) {
    ## Cellassign
    sobj <- cell.annotation.cellassign(sobj = sobj, markers = cellassign.markers, GPU = FALSE, by.cluster = TRUE)
    
    # Check markers inside the annotated cluster  
    # markers_SCT <- Seurat::DotPlot(object = sobj, features = cellassign.markers, assay = 'SCT', group.by = 'cellassign') + 
    #   ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
    markers_RNA <- Seurat::DotPlot(object = sobj, features = unique(cellassign.markers), assay = 'RNA', group.by = 'cellassign') + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
    
    # Plot the annotation on UMAP
    cellassign_plot <- Seurat::DimPlot(sobj,group.by = "cellassign",reduction = umap.name, pt.size = solo.pt.size) + 
      ggplot2::ggtitle(paste0("Cellassign cell type prediction for ",ncol(sobj), "cells")) +
      ggplot2::theme(legend.position = "bottom") + 
      Seurat::DarkTheme()
    cellassign_cluster_plot <- Seurat::DimPlot(sobj,group.by = "cellassign_cluster",reduction = umap.name, pt.size = solo.pt.size) + 
      ggplot2::ggtitle(paste0("Cellassign cell type prediction for ",nlevels(sobj@meta.data[['cellassign_cluster']]), "cluster(s)")) +
      ggplot2::theme(legend.position = "bottom") + 
      Seurat::DarkTheme()
    
    ggplot2::ggsave(filename = paste0(sample.name, '_', assay, '_uMAPs_', "cellasign", '.png'), path = cellassign.dir, plot= cellassign_plot, width = 8, height = 9)
    ggplot2::ggsave(filename = paste0(sample.name, '_', assay, '_uMAPs_', "cellasign_clust", '.png'), path = cellassign.dir, plot= cellassign_cluster_plot, width = 8, height = 9)
    # ggplot2::ggsave(filename = paste0(sample.name, '_', assay, '_uMAPs_', "markers_SCT_check", '.png'), path = cellassign.dir, plot= markers_SCT, width = 4, height = 4)
    ggplot2::ggsave(filename = paste0(sample.name, '_', assay, '_uMAPs_', "markers_RNA_check", '.png'), path = cellassign.dir, plot= markers_RNA, width = 10, height = 5)
    sobj@misc$params$cellassign.markers <- cellassign.markers
  }
  return(sobj)
}

## Find markers
find.markers.quick <- function(sobj = NULL, assay = 'SCT', ident = NULL, slot, test.use = 'wilcox', min.pct = .75, logfc.threshold = .5, only.pos = TRUE, adjp.p.max = 1E-02, topn = 10, heatmap.cols = c("gold", "blue"), out.dir = NULL) {
  
  if (is.null(out.dir)) stop('No output dir provided !')
  if (!dir.exists(out.dir)) stop('Output directory does not exist !')
  if (is.null(sobj)) stop('No Seurat object provided !')
  if (!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist !'))
  if (!"seed" %in% names(sobj@misc$params)) stop("No seed found in Seurat object !")
  if (!is.null(ident)) {
    if (!ident %in% names(sobj@meta.data)) stop(paste0("Could not find '", ident, "' in @meta.data !"))
  } else stop("No ident name provided !")
  
  my.seed <- sobj@misc$params$seed
    
  sample.name <- Seurat::Project(sobj)
  
  suppressPackageStartupMessages(require(UpSetR))
  suppressPackageStartupMessages(require(dplyr))
  
  if (sum(dim(sobj@assays[[assay]]@scale.data)) == 0) {
    sobj <- Seurat::ScaleData(object = sobj, assay = assay, features = rownames(sobj@assays[[assay]]@counts))
  }
  
  ## Keeping track of original ident (to restore it after our computations)
  ori.ident <- Seurat::Idents(sobj)
  Seurat::Idents(sobj) <- ident
  
  ## Computing differentials
  fmark <- Seurat::FindAllMarkers(sobj, assay = assay, test.use = test.use, min.pct = min.pct, logfc.threshold = logfc.threshold, only.pos = only.pos, random.seed = my.seed)
  
  ## Filtering markers
  fmark <- fmark[fmark$p_val_adj < adjp.p.max, ]
  fmark <- fmark[order(fmark$cluster, fmark$p_val_adj),]
  mytop <- fmark %>% group_by(cluster) %>% top_n(n = topn, wt = avg_logFC)
  
  ## Creating output dir
  fmark.dir <- paste0(out.dir, '/found_markers_', ident, '/')
  dir.create(fmark.dir, recursive = TRUE, showWarnings = FALSE)
  
  ## Heatmap
  suppressPackageStartupMessages(require(ggplot2))
  png(paste0(fmark.dir, '/', sample.name, '_findmarkers_top', topn, '_heatmap.png'), width = 1600, height = 1000)
  print(Seurat::DoHeatmap(sobj, features = mytop$gene, angle = 0, hjust = .5, assay = assay) + ggplot2::scale_fill_gradientn(colors = heatmap.cols))
  dev.off()
  
  ## Upset plots
  g.list <- sapply(sort(unique(fmark$cluster)), function(k) { return(fmark$gene[fmark$cluster == k]) }, simplify = FALSE)
  names(g.list) <- sort(unique(fmark$cluster))
  g.top.list <- sapply(sort(unique(mytop$cluster)), function(k) { return(mytop$gene[mytop$cluster == k]) }, simplify = FALSE)
  names(g.top.list) <- sort(unique(mytop$cluster))
  png(paste0(fmark.dir, '/', sample.name, '_findmarkers_upset_all.png'), width = 1600, height = 1000)
  up1 <- upset(fromList(g.list), sets = rev(names(g.list)), nsets = length(g.list), text.scale = 3, keep.order = TRUE)
  print(up1)
  grid::grid.text("All markers",x = 0.65, y=0.95, gp=grid::gpar(fontsize=30))
  dev.off()
  png(paste0(fmark.dir, '/', sample.name, '_findmarkers_upset_top', topn, '.png'), width = 1600, height = 1000)
  up2 <- upset(fromList(g.top.list), sets = rev(names(g.top.list)), nsets = length(g.top.list), text.scale = 3, keep.order = TRUE)
  print(up2)
  grid::grid.text(paste0("Top ", topn, " markers"),x = 0.65, y=0.95, gp=grid::gpar(fontsize=30))
  dev.off()
  
  ## Violinplots
  for(k in unique(fmark$cluster)) {
    mytop.k <- mytop[mytop$cluster == k,]
    if (nrow(mytop.k) > 0) {
      # Seurat::VlnPlot(sobj, features = mytop.k$gene, ncol = 5)
      fm.list <- sapply(seq_len(nrow(mytop.k)), function(g) { Seurat::VlnPlot(sobj, features = mytop.k$gene[g]) + Seurat::NoLegend()+ ggplot2::ggtitle(paste0(mytop.k$gene[g], ' ; logFC = ', format(mytop.k$avg_logFC[g], digits= 3), ' ; adj.p = ', format(mytop.k$p_val_adj[g], digits = 2, scientific = TRUE))) }, simplify = FALSE)
      png(paste0(fmark.dir, '/', sample.name, '_findmarkers_top', topn, '_cluster', k, '_vln.png'), width = 1800, height = 1000)
      # print(Seurat::CombinePlots(plots = fm.list, ncol = 5))
      print(patchwork::wrap_plots(fm.list) + patchwork::plot_layout(ncol = 3))
      dev.off()
    }
  }
  
  ## Restoring original ident
  Seurat::Idents(sobj) <- ori.ident
  
  
  ## Save table
  write.table(fmark, file = paste0(fmark.dir, '/', sample.name, '_', ident, '_findmarkers_all.txt'), sep = "\t", row.names = FALSE, quote = FALSE)
  
  sobj@assays[[assay]]@scale.data <- matrix(nrow = 0, ncol = 0)
  
  sobj@misc$find.markers.quick <- fmark
  sobj@misc$params$find.markers.quick$test.use = test.use
  sobj@misc$params$find.markers.quick$min.pct = min.pct
  sobj@misc$params$find.markers.quick$logfc.threshold = logfc.threshold
  sobj@misc$params$find.markers.quick$only.pos = TRUE
  sobj@misc$params$find.markers.quick$topn = 10
  sobj@misc$Rsession <- sessionInfo()

  return(sobj) 
}

## Technical uMAPs
technical.plot <- function(sobj = NULL, assay = 'RNA', umap.name = NULL, ident = NULL, out.dir = NULL, multi.pt.size = 1) {
  
  if (is.null(out.dir)) stop('No output dir provided !')
  if (!dir.exists(out.dir)) stop('Output directory does not exist !')
  if (is.null(sobj)) stop('No Seurat object provided !')
  if (is.null(umap.name)) stop('No umap reduction name provided !')
  
  tech.dir <- paste0(out.dir, '/technical/')
  dir.create(tech.dir, recursive = TRUE, showWarnings = TRUE)
  
  require(patchwork)
  
  sample.name <- Seurat::Project(sobj)
  
  plot.num.add <- if(!is.null(ident)) 1 else 0
  plot.pix <- 600
  
  ## CLUSTERS
  if (!is.null(ident)) {
    # ori.ident <- Seurat::Idents(sobj)
    Seurat::Idents(sobj) <- ident
    upCLUST <- Seurat::LabelClusters(plot = Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size) + ggplot2::ggtitle("Louvain clusters (Seurat)") + Seurat::DarkTheme(), id = "ident", size = multi.pt.size*3, repel = FALSE, color = "white", fontface = "bold")
    # Seurat::Idents(sobj) <- ori.ident
  }
  
  ## METRICS
  metrics.plotlist <- list(
    'percent.mb' = Seurat::FeaturePlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, features = "percent.mt", max.cutoff = .1, cols = gradient.cols, order = TRUE) + Seurat::DarkTheme(),
    'percent.rb' = Seurat::FeaturePlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, features = "percent.rb", max.cutoff = .25, cols = gradient.cols, order = TRUE) + Seurat::DarkTheme(),
    'log_nCount_RNA' = Seurat::FeaturePlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, features = "log_nCount_RNA", cols = gradient.cols, order = TRUE) + Seurat::DarkTheme(),
    'nFeature_RNA' = Seurat::FeaturePlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, features = "nFeature_RNA", cols = gradient.cols, order = TRUE) + Seurat::DarkTheme()
  )
  
  for (p in seq_along(metrics.plotlist)) {
    png(paste0(tech.dir, '/', sample.name, '_technical_SINGLE_METRICS_', names(metrics.plotlist)[p], '_uMAP.png'), width = plot.pix, height = plot.pix)
    print(metrics.plotlist[[p]])
    dev.off()
  }
  grid.xy <- grid.scalers(length(metrics.plotlist) + plot.num.add)
  png(paste0(tech.dir, '/', sample.name, '_technical_MULTI_METRICS_uMAPs.png'), width = grid.xy[1]*plot.pix, height = grid.xy[2]*plot.pix)
  if(!is.null(ident)) print(upCLUST + metrics.plotlist) else print(metrics.plotlist) 
  dev.off()
  
  ## CELL CYCLE (cyclone)
  cycle.plotlist <- list()  
  if ('Phase' %in% colnames(sobj@meta.data)) {
    cycle.plotlist <- list(
      'cyclone_Phase' = Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, group.by = "Phase")  + ggplot2::ggtitle("Cell Phase (cyclone)") + Seurat::DarkTheme(),
      'cyclone_S' = Seurat::FeaturePlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, features = "S.Score", cols = gradient.cols, order = TRUE) + Seurat::DarkTheme(),
      'cyclone_G2M' = Seurat::FeaturePlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, features = "G2M.Score", cols = gradient.cols, order = TRUE) + Seurat::DarkTheme(),
      'cyclone_SmG2M' = Seurat::FeaturePlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, features = "SmG2M.Score", cols = gradient.cols, order = TRUE) + Seurat::DarkTheme(),
      'cyclone_G1' = Seurat::FeaturePlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, features = "G1.Score", cols = gradient.cols, order = TRUE) + Seurat::DarkTheme()
    )
    for (p in seq_along(cycle.plotlist)) {
      png(paste0(tech.dir, '/', sample.name, '_technical_SINGLE_CYCLE_', names(cycle.plotlist)[p], '_uMAP.png'), width = 1000, height = 1000)
      print(cycle.plotlist[[p]])
      dev.off()
    }
    grid.xy <- grid.scalers(length(cycle.plotlist) + plot.num.add)
    png(paste0(tech.dir, '/', sample.name, '_technical_MULTI_CYCLE_cyclone_uMAPs.png'), width = grid.xy[1]*plot.pix, height = grid.xy[2]*plot.pix)
    if(!is.null(ident)) print(upCLUST + cycle.plotlist) else print(cycle.plotlist) 
    dev.off()
  }
  
  ## CELL CYCLE (Seurat)
  s_cycle.plotlist <- list()  
  if ('Seurat.Phase' %in% colnames(sobj@meta.data)) {
    s_cycle.plotlist <- list(
      'Seurat_Phase' = Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, group.by = "Seurat.Phase")  + ggplot2::ggtitle("Cell Phase (Seurat)") + Seurat::DarkTheme(),
      'Seurat_S' = Seurat::FeaturePlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, features = "Seurat.S.Score", cols = gradient.cols, order = TRUE) + Seurat::DarkTheme(),
      'Seurat_G2M' = Seurat::FeaturePlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, features = "Seurat.G2M.Score", cols = gradient.cols, order = TRUE) + Seurat::DarkTheme(),
      'Seurat_SmG2M' = Seurat::FeaturePlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, features = "Seurat.SmG2M.Score", cols = gradient.cols, order = TRUE) + Seurat::DarkTheme()
    )
    for (p in seq_along(s_cycle.plotlist)) {
      png(paste0(tech.dir, '/', sample.name, '_technical_SINGLE_CYCLE_', names(s_cycle.plotlist)[p], '_uMAP.png'), width = 1000, height = 1000)
      print(s_cycle.plotlist[[p]])
      dev.off()
    }
    grid.xy <- grid.scalers(length(s_cycle.plotlist) + plot.num.add)
    png(paste0(tech.dir, '/', sample.name, '_technical_MULTI_CYCLE_Seurat_uMAPs.png'), width = grid.xy[1]*plot.pix, height = grid.xy[2]*plot.pix)
    if(!is.null(ident)) print(upCLUST + s_cycle.plotlist) else print(s_cycle.plotlist) 
    dev.off()
  }
  
  ## CELL CYCLE (all)
  if (length(cycle.plotlist) > 0 & length(s_cycle.plotlist) > 0) {
    grid.xy <- grid.scalers(length(cycle.plotlist) + length(s_cycle.plotlist) + plot.num.add)
    png(paste0(tech.dir, '/', sample.name, '_technical_MULTI_CYCLE_ALL_uMAPs.png'), width = grid.xy[1]*plot.pix, height = grid.xy[2]*plot.pix)
    if(!is.null(ident)) print(upCLUST + cycle.plotlist + s_cycle.plotlist) else print(cycle.plotlist + s_cycle.plotlist) 
    dev.off()
  }
  
  
  ## DOUBLETS
  doublets.plotlist <- list()
  if("scDblFinder.class" %in% colnames(sobj@meta.data)) doublets.plotlist[['ScDblFinder']] <- Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, group.by = "scDblFinder.class") + ggplot2::ggtitle("Cell doublets (scDblFinder)") + Seurat::DarkTheme()
  if("hybrid_score.class" %in% colnames(sobj@meta.data)) doublets.plotlist[['scds_hybrid']] <- Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, group.by = "hybrid_score.class") + ggplot2::ggtitle("Cell doublets (scds-hybrid)") + Seurat::DarkTheme()
  if("doublets_consensus.class" %in% colnames(sobj@meta.data)) doublets.plotlist[['Doublets_union']] <- Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, group.by = "doublets_consensus.class") + ggplot2::ggtitle("Cell doublets (union)") + Seurat::DarkTheme()
  if("log_scran.doubletscore" %in% colnames(sobj@meta.data)) doublets.plotlist[['Doublets_scran']] <- Seurat::FeaturePlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, features = "log_scran.doubletscore", cols = gradient.cols, order = TRUE)  + ggplot2::ggtitle("Cell doublets (scran (log))") + Seurat::DarkTheme()
  
  if (length(doublets.plotlist) > 0) {
    for (p in seq_along(doublets.plotlist)) {
      png(paste0(tech.dir, '/', sample.name, '_technical_SINGLE_', names(doublets.plotlist)[p], '_uMAP.png'), width = 1000, height = 1000)
      print(doublets.plotlist[[p]])
      dev.off()
    }
    grid.xy <- grid.scalers(length(doublets.plotlist) + plot.num.add)
    png(paste0(tech.dir, '/', sample.name, '_technical_MULTI_DOUBLETS_uMAPs.png'), width = grid.xy[1]*plot.pix, height = grid.xy[2]*plot.pix)
    if(!is.null(ident)) print(upCLUST + doublets.plotlist) else print(doublets.plotlist) 
    dev.off()
  }

  ## MULTI:ALL
  grid.xy <- grid.scalers(length(metrics.plotlist) + length(doublets.plotlist) + length(cycle.plotlist) + length(s_cycle.plotlist) + plot.num.add)
  png(paste0(tech.dir, '/', sample.name, '_technical_MULTI_ALL_uMAPs.png'), width = grid.xy[1]*plot.pix, height = grid.xy[2]*plot.pix)
  if(!is.null(ident)) print(upCLUST + metrics.plotlist + doublets.plotlist + cycle.plotlist + s_cycle.plotlist) else print(metrics.plotlist + doublets.plotlist + cycle.plotlist + s_cycle.plotlist) 
  dev.off()
}

## Markers
markers.umap.plot <- function(sobj = NULL, markers = NULL, assay = 'RNA', umap.name = NULL, ident = NULL, out.dir = NULL, dimplot.cols = c("gold", "blue"), multi.pt.size = 1) {
  if (is.null(out.dir)) stop('No output dir provided !')
  if (!dir.exists(out.dir)) stop('Output directory does not exist !')
  if (is.null(sobj)) stop('No Seurat object provided !')
  if (is.null(umap.name)) stop('No UAMP reduction name provided !')
  if (!(umap.name %in% names(sobj@reductions))) stop(paste0('UMAP name "', umap.name, '" does not exist !'))
  if (!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist !'))
  if (is.null(markers)) stop('No marker name provided !')
  
  sample.name <- Seurat::Project(sobj)
  
  markers <- markers[markers %in% rownames(sobj@assays[[assay]]@data)]
  
  mark.dir <- paste0(out.dir, '/markers/')
  dir.create(mark.dir, recursive = TRUE, showWarnings = FALSE)
  
  require(patchwork)
  
  plot.num.add <- if(!is.null(ident)) 1 else 0
  plot.pix <- 600
  
  ## CLUSTERS
  if (!is.null(ident)) {
    # ori.ident <- Seurat::Idents(sobj)
    Seurat::Idents(sobj) <- ident
    upCLUST <- Seurat::LabelClusters(plot = Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size) + ggplot2::ggtitle("Louvain clusters (Seurat)") + Seurat::DarkTheme(), id = "ident", size = multi.pt.size*3, repel = FALSE, color = "white", fontface = "bold")
    # Seurat::Idents(sobj) <- ori.ident
  }
  
  marklist <- list()
  marknames <- unique(names(markers))
  if (length(marknames) > 0) {
    for (mn in marknames) {
      mini.markers <- markers[names(markers) == mn]
      mn.plotlist <- sapply(mini.markers, function(x) {Seurat::FeaturePlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, features = x, cols = dimplot.cols) + Seurat::DarkTheme() }, simplify = FALSE)
      grid.xy <- grid.scalers(length(mn.plotlist) + plot.num.add)
      png(paste0(mark.dir, '/', sample.name, '_markers_TYPE_', gsub(pattern = " ", replacement = "_", mn), '_uMAPs.png'), width = grid.xy[1]*plot.pix, height = grid.xy[2]*plot.pix)
      if(!is.null(ident)) print(upCLUST + mn.plotlist) else print(mn.plotlist) 
      dev.off()
      
      marklist <- c(marklist, list(mn.plotlist))
    }
  }
  marklist <- unlist(marklist, recursive = FALSE)
  
  ## SINGLES
  for(x in seq_along(marklist)) {
    png(paste0(mark.dir, '/', sample.name, '_markers_SINGLE_', markers[x], '_uMAP.png'), width = 1000, height = 1000)
    print(marklist[[x]])
    dev.off()
  }
  ## ALL
  if (length(markers) <= 15) {
    grid.xy <- grid.scalers(length(marklist) + plot.num.add)
    png(paste0(mark.dir, '/', sample.name, '_markers_ALL_uMAPs.png'), width = grid.xy[1]*plot.pix, height = grid.xy[2]*plot.pix)
    if(!is.null(ident)) print(upCLUST + marklist) else print(marklist) 
    dev.off()
  }
  sobj@misc$markers.in <- markers
  return(sobj)
}

## Control genes
ctrl.umap.plot <- function(sobj = NULL, ctrl.genes = NULL, assay = 'RNA', umap.name = NULL, ident = NULL, out.dir = NULL, dimplot.cols = c("gold", "blue"), multi.pt.size = 2) {
  if (is.null(out.dir)) stop('No output dir provided !')
  if (!dir.exists(out.dir)) stop('Output directory does not exist !')
  if (is.null(sobj)) stop('No Seurat object provided !')
  if (is.null(umap.name)) stop('No UAMP reduction name provided !')
  if (!(umap.name %in% names(sobj@reductions))) stop(paste0('UMAP name "', umap.name, '" does not exist !'))
  if (!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist !'))
    
  sample.name <- Seurat::Project(sobj)
  
  
  ctrl.genes <- unique(ctrl.genes[ctrl.genes %in% rownames(sobj@assays[[assay]]@data)])
  sobj@misc$ctrl.genes <- ctrl.genes
  
  plot.num.add <- if(!is.null(ident)) 1 else 0
  plot.pix <- 600
  
  ## CLUSTERS
  if (!is.null(ident)) {
    # ori.ident <- Seurat::Idents(sobj)
    Seurat::Idents(sobj) <- ident
    upCLUST <- Seurat::LabelClusters(plot = Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size) + ggplot2::ggtitle("Louvain clusters (Seurat)") + Seurat::DarkTheme(), id = "ident", size = multi.pt.size*3, repel = FALSE, color = "white", fontface = "bold")
    # Seurat::Idents(sobj) <- ori.ident
  }
  
  ctrl.dir <- paste0(out.dir, '/control_genes/')
  dir.create(ctrl.dir, recursive = TRUE, showWarnings = FALSE)
  
  require(patchwork)
  
  for (x in ctrl.genes) {
    exp.umap.counts <- Seurat::FeaturePlot(object = sobj, reduction = umap.name, slot = 'counts', pt.size = multi.pt.size, features = x, cols = dimplot.cols) + ggplot2::ggtitle(paste0(x, ' expression level (counts)')) + Seurat::DarkTheme()
    exp.umap.log <- Seurat::FeaturePlot(object = sobj, reduction = umap.name, slot = 'data', pt.size = multi.pt.size, features = x, cols = dimplot.cols) + ggplot2::ggtitle(paste0(x, ' expression level (logcounts)')) + Seurat::DarkTheme()
    ctrl.list <- list(exp.umap.counts, exp.umap.log)
    if(paste0('ctrl_', assay, '_', x) %in% colnames (sobj@meta.data)) {
      tag.umap <- Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, group.by = paste0('ctrl_', assay, '_', x)) + ggplot2::ggtitle(paste0(x, ' is expressed')) + Seurat::DarkTheme()
      ctrl.list <- append(ctrl.list, list(tag.umap))
    }
    
    grid.xy <- grid.scalers(length(ctrl.list) + plot.num.add)
    png(paste0(ctrl.dir, '/', sample.name, '_control.genes_', x, '_uMAP.png'), width = grid.xy[1]*plot.pix, height = grid.xy[2]*plot.pix)
    if(!is.null(ident)) print(upCLUST + ctrl.list) else print(ctrl.list) 
    dev.off()
  }
  sobj@misc$ctrl.genes <- ctrl.genes
  return(sobj)
}

## Build annotated Cerebro binary
### Additional parameters (...) correspond to the cerebroApp::getMarkerGenes() function, itself forwarding parameters to the Seurat::FindAllMarkers() function
### Those recommended parameters are : only_pos = TRUE, min_pct = .75, thresh_logFC = .5, thresh_p_val = 1E-02, test = 'wilcox'
seurat2cerebro <- function(sobj = NULL, assay = 'RNA', ident = NULL, species = 'hg', gmt.file = 'msigdb.v7.1.symbols.gmt', out.dir = getwd(), nthreads = 1, ...) {
  if (is.null(sobj)) stop('A seurat object is required !')
  if (is.null(ident)) stop('An ident is required (as default clustering) !')
  if (!ident %in% colnames(sobj@meta.data)) stop('ident not found in Seurat object !')
  if (!species %in% c('hg', 'mm', 'rn'))
  sample.name <- Seurat::Project(sobj)
  sobj$cluster <- sobj[[ident]]
  Idents(sobj) <- sobj$cluster
  sobj$nUMI <- sobj[[paste0('nCount_', assay)]]
  sobj$nGene <- sobj[[paste0('nFeature_', assay)]]
  sobj <- cerebroApp::getMostExpressedGenes(object = sobj, assay = assay)
  sobj <- cerebroApp::getMarkerGenes(object = sobj, assay = assay, organism = species, ...)
  sobj <- cerebroApp::getEnrichedPathways(object = sobj)
  sobj <- cerebroApp::performGeneSetEnrichmentAnalysis(object = sobj, assay = assay, GMT_file = gmt.file, parallel.sz = nthreads)
  cerebroApp::exportFromSeurat(object = sobj, assay = assay, experiment_name = paste0(sample.name, '_', assay), organism = species, column_cell_cycle_seurat = 'Seurat.Phase', column_cell_cycle_cyclone = 'Phase', file = paste0(out.dir, '/', paste(c(sample.name, assay, ident), collapse = '_'), '.crb'))
  return(sobj)
}

## Computing correlations between two assays with corresponding features (ex: RNA and ADT)
### assay features must be in the same order !
feature.cor <- function(sobj = NULL, assay1 = 'RNA', assay2 = 'ADT', assay1.features = NULL, assay2.features = NULL, slot = 'data', cor.method = 'spearman', zero.filter = TRUE) {
  
  cell.idx.list <- sapply(
    seq_along(gene.names),
    function(k) {
      if (zero.filter) zerocells.get(sobj = sobj, assay1 = assay1, assay2 = assay2, assay1.feature = assay1.features[k], assay2.feature = assay2.features[k], slot = slot) else rep(TRUE, ncol(sobj@assays[[assay1]]))
    },
    simplify = FALSE
  )
    
  corvec <- vapply(
    seq_along(gene.names),
    function(k) {
      cell.idx <- cell.idx.list[[k]]
      # if(!any(cell.idx)) return(NA)
      if(length(which(cell.idx)) < 2) return(NA)
      return(
        cor.test(
          slot(sobj@assays[[assay1]], slot)[rownames(slot(sobj@assays[[assay1]], slot)) == assay1.features[k],cell.idx],
          slot(sobj@assays[[assay2]], slot)[rownames(slot(sobj@assays[[assay2]], slot)) == assay2.features[k],cell.idx],
          method = cor.method
        )$estimate
      )
    },
  .1)
  out.df <- if (zero.filter) {
    data.frame(vapply(cell.idx.list, function(x) { length(which(x)) }, 1L),  corvec, stringsAsFactors = FALSE)
  } else {
    data.frame(corvec, stringsAsFactors = FALSE)
  }
  colnames(out.df) <- if (zero.filter) {
    c(paste(c(slot, 'non0'), collapse = '_'), paste(c('cor', slot, cor.method, '0filt'), collapse = '_'))
  } else {
    c(paste(c('cor', slot, cor.method), collapse = '_'))
  }
  
  return(out.df)
  # if (zero.filter) {
  #   out.df <- data.frame(vapply(cell.idx.list, function(x) { length(which(x)) }, 1L),  corvec, stringsAsFactors = FALSE)
  #   colnames(out.df) <- c(paste(c('cor', slot, 'non0'), collapse = '_'), paste(c('cor', slot, cor.method, '0filt'), collapse = '_'))
  # } else {
  #   out.df <-data.frame(corvec, stringsAsFactors = FALSE)
  #   colnames(out.df) <- c(paste(c('cor', slot, cor.method, '0filt'), collapse = '_'))
  # }
}

zerocells.get <- function(sobj = NULL, assay1 = 'RNA', assay2 = 'ADT', assay1.feature = NULL, assay2.feature = NULL, slot = 'data') {
  cell.idx1 <- slot(sobj@assays[[assay1]], slot)[rownames(slot(sobj@assays[[assay1]], slot)) == assay1.feature,] > 0
  cell.idx2 <- slot(sobj@assays[[assay2]], slot)[rownames(slot(sobj@assays[[assay2]], slot)) == assay2.feature,] > 0
  return(cell.idx1 & cell.idx2)
}

## Get the list of assays of a Seurat object
get.assays <- function(sobj = NULL) {
  if (is.null(sobj)) stop('No Seurat object provided !')
  print(Seurat::Assays(sobj))
}

## Get the list of reductions of a Seurat object
get.reductions <- function(sobj = NULL) {
  if (is.null(sobj)) stop('No Seurat object provided !')
  print(Seurat::Reductions(sobj))
}

## Get the list of idents that can be used out of a Seurat object
get.idents <- function(sobj = NULL, pattern = '_res') {
  if (is.null(sobj)) stop('No Seurat object provided !')
  idents.list <- grep(pattern = pattern, x = colnames(sobj@meta.data))
  for (id in idents.list) {
    message(paste0("\n", colnames(sobj@meta.data)[id]))
    print(table(sobj@meta.data[,id]))
  }
}

## Get the list of graphs of a Seurat object
get.graphs <- function(sobj = NULL) {
  if (is.null(sobj)) stop('No Seurat object provided !')
  # print(grep("_res", x = colnames(sobj@meta.data), value = T))
  print(names(sobj@graphs))
}


## Function to get the multiplot scalers
grid.scalers <- function(n = 1) {
  x <- ceiling(sqrt(n))
  y <- ceiling(n / x)
  # print(paste0("X : ", x))
  # print(paste0("Y : ", y))
  return(c(x, y))
}

load_ressources <- function(path, species){
  supported_species <- fs::path(path, "GENELISTS") %>% fs::dir_ls(type = "directory") %>% fs::path_file()
  assertthat::assert_that(is_in(species, supported_species))
  
  if (species == "homo_sapiens") {
    tp <- list(species.rdx = 'hg',
    mt.genes.file = fs::path(path, "GENELISTS/homo_sapiens/homo_sapiens_mito_symbols_20191001.rds"),
    crb.genes.file = fs::path(path, "GENELISTS/homo_sapiens/homo_sapiens_cribo_symbols_20191015.rds"),
    cc.genes.file = fs::path(path, "GENELISTS/homo_sapiens/homo_sapiens_Seurat_cc.genes_20191031.rds"),
    cc.pairs.file = fs::path(path, "GENELISTS/homo_sapiens/homo_sapiens_cyclone_pairs_symbols_20191001.rds"),
    str.genes.file = fs::path(path, "GENELISTS/homo_sapiens/homo_sapiens_stress_symbols_20200224.rds"),
    singler.setnames = c("HumanPrimaryCellAtlasData", "NovershternHematopoieticData", "DatabaseImmuneCellExpressionData", "MonacoImmuneData"),
    clustifyr.setnames = c("pbmc_avg", "hema_microarray_matrix", "gtex_bulk_matrix"),
    gmt.file = fs::path(path, '/RESOURCES/DATABASE/MSIGDB/7.1/msigdb_v7.1_GMTs/msigdb.v7.1.symbols.gmt')) %>%
      list2env(envir = .GlobalEnv)
  } else if (species == "mus_musculus") {
    tp <- list(species.rdx = 'mm',
    mt.genes.file = fs::path(path, "GENELISTS/mus_musculus/mus_musculus_mito_symbols_20191015.rds"),
    crb.genes.file = fs::path(path, "GENELISTS/mus_musculus/mus_musculus_cribo_symbols_20191015.rds"),
    cc.genes.file = fs::path(path, "GENELISTS/mus_musculus/mus_musculus_Seurat_cc.genes_20191031.rds"),
    cc.pairs.file = fs::path(path, "GENELISTS/mus_musculus/mus_musculus_cyclone_pairs_symbols_20191015.rds"),
    str.genes.file = fs::path(path, "GENELISTS/mus_musculus/mus_musculus_stress_symbols_20200224.rds"),
    singler.setnames = c("MouseRNAseqData", "ImmGenData"),
    clustifyr.setnames = c("ref_MCA", "ref_tabula_muris_drop", "ref_tabula_muris_facs", "ref_moca_main", "ref_immgen", "ref_mouse.rnaseq"),
    gmt.file = fs::path(path, 'msigdb/v7.1/msigdb.v7.1.symbols.gmt')) %>%
      list2env(envir = .GlobalEnv)
  }
}

read_markers <- function(path, sheet = "MARKERS"){
  openxlsx::read.xlsx(path, sheet = sheet, startRow = 1, fillMergedCells = TRUE, colNames = TRUE) %>%
    tibble::rowid_to_column() %>%
    tidyr::pivot_longer(cols = -rowid, values_drop_na = TRUE) %>%
    dplyr::select(-rowid) %>%
    dplyr::arrange(name) %>%
    tibble::deframe()
}

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


## DATASETS
# rootdir <- "/mounts/cigogne/data/job/SingleCell/scRNAseq/TESTS/KALLISTO_BUSTOOLS"
# rootdir <- "/media/job/dcc7943e-b81b-4c6e-9b2f-9159b58c7824/job/SINGLECELL/"
# rootdir <- "/mounts/flamingo/beegfs/scratch/SCRNASEQ/"
## Sourcing functions
source("Clustering//bustools2seurat_preproc_functions.R")
source("utils.R")
source("Clustering/cellassign.R")


rootdir <- "/mnt/Data/Project/"
rundir <- "/mnt/Data"
run.name <- "Project"

project.name <- "NF1-MPNST_Human"
sample.name <- "2020_12"

species <- "homo_sapiens"
## The range of mitochondrial gene expression percentage to keep cell in.
pcmito.range <- c(0, .1)
## The range of ribosomal gene expression percentage to keep cell in.
pcribo.range <- c(0, 1)
## Minimum features (genes) to keep a cell
min.features <- 200
## Minimum counts (UMIs) to keep a cell
min.counts <- 1000
## Minimum cells with expression > 0 to keep a gene
min.cells <- 5
## Minimal cluster size for overclustering steps (see scran/scater help)
min.clust.size <- 50
markers <- read_markers("/mnt/Data/Project/test.markers.xlsx", "MARKERS")
markers.cellassign <- read_markers("/mnt/Data/Project/test.markers.xlsx", "CELLASSIGN")
## PARAMETERS
emptydrops.fdr <- 1E-03  ## FDR to select unique droplets
cfr.minscore <- 0.35 ## Minimum correlation score for clustifyr to consider
sr.minscore <- 0.25 ## Minimum correlation score for SingleR to consider
droplets.limit <- 1E+05 ## Max number of droplets to activate empty droplets filtering
red.dims <- 100L
my.seed <- 1337L
nthreads <- 6L
solo.pt.size <- 3L
multi.pt.size <- 2L
gradient.cols <- c("gold", "blue")
gradient.cols3 <- c("gold", "white", "blue")
ressources_path <- "/mnt/Data/RESSOURCES/"
ctrl.genes <- "dtTomato"
eval.markers <- c("Gadph", ctrl.genes)
## Gene lists loading
load_ressources(ressources_path, species)



## SAMPLE/DIR
data.path <- paste0(rundir, "/", run.name, "/", project.name, "/", sample.name, "/")
dir.exists(data.path)


## RUN! (EXAMPLES)
#######

## A. UNFILTERED ROLL
##-------------------

### PARAMS
assay <- 'RNA'
norm.method <- 'LogNormalize' ## or 'SCTransform'
reduction <- 'pca'
keep.dims <- 15 ## forced
keep.res <- .8 ## forced

### Creating parallel instance
cl <- create.parallel.instance(nthreads = nthreads)

### Loading raw count matrix
sobj <- load.sc.data(data.path = data.path, sample.name = sample.name, droplets.limit = droplets.limit, emptydrops.fdr = emptydrops.fdr, BPPARAM = cl, my.seed = my.seed)

### Computing basic metrics
sobj <- QC.metrics(sobj = sobj, mt.genes.file = mt.genes.file, crb.genes.file = crb.genes.file, str.genes.file = str.genes.file, pcmito.range = pcmito.range, pcribo.range = pcribo.range, min.features = min.features, min.counts = min.counts, BPPARAM = cl)

### Tag control genes
sobj <- tag.ctrl.genes(sobj = sobj, ctrl.genes = ctrl.genes, ctrl.min.counts = 3, assay = "RNA")

### Building output directiry
out.dir <- paste0(data.path, '/UNFILTERED/')
dir.create(path = out.dir, recursive = TRUE, showWarnings = TRUE)

### Saving non-normalized object
save(sobj, file = paste0(out.dir, '/', sample.name, '_UNFILTERED_NON-NORMALIZED.rda'), compress = "bzip2")

### Basic normalization
sobj <- sc.normalization(sobj = sobj, assay = assay, method = norm.method, features.n = 3000)

### Building normalized output dir
norm.dir <- paste(out.dir, norm.method, sep = '/')
dir.create(path = norm.dir, recursive = TRUE, showWarnings = TRUE)

### Saving normalized object
save(sobj, file = paste0(norm.dir, '/', sample.name, '_UNFILTERED_', norm.method, '.rda'), compress = "bzip2")

### Dim reduction
sobj <- dimensions.reduction(sobj = sobj, reduction = reduction, assay = assay, max.dims = keep.dims)

### Building clustered output dir
clust.dir <- paste(norm.dir, paste(c(reduction, keep.dims, keep.res), collapse = '_'), sep = '/')
dir.create(path = clust.dir, recursive = TRUE, showWarnings = TRUE)

## Clustering + uMAP
sobj <- louvain.cluster(sobj = sobj, reduction = reduction, assay = assay, max.dim = keep.dims, resolution = keep.res, out.dir = clust.dir, solo.pt.size = solo.pt.size)

## Technical plots
technical.plot(sobj = sobj, assay = assay, umap.name = paste(c(assay, reduction, keep.dims, 'umap'), collapse = '_'), ident = paste0(paste(c(assay, reduction, keep.dims), collapse = '.'), '_res.', keep.res), out.dir = clust.dir, multi.pt.size = multi.pt.size)

### Saving reduced object
save(sobj, file = paste0(clust.dir, '/', paste(c(sample.name, 'UNFILTERED', norm.method, reduction, keep.dims, keep.res), collapse = '_'), '.rda'), compress = "bzip2")


## B. FILTERED ROLL, DOUBLETS KEPT
##--------------------------------

### PARAMS
assay <- 'RNA'
norm.method <- 'SCTransform' ## or 'SCTransform'
reduction <- 'pca'
keep.dims <- 15 ## forced
keep.res <- .8 ## forced
vtr <- NULL ## No covariate to regress

expected.rda <- paste0(
  paste(
    paste0(data.path, '/UNFILTERED/'),
    norm.method,
    sep = '/'
  ),
  '/',
  paste(c(sample.name, 'UNFILTERED', norm.method), collapse = '_'),
  '.rda'
)

if (file.exists(expected.rda)) {
  load(expected.rda)
} else {
  ### Creating parallel instance
  cl <- create.parallel.instance(nthreads = nthreads)
  
  ### Loading raw count matrix
  sobj <- load.sc.data(data.path = data.path, sample.name = sample.name, droplets.limit = droplets.limit, emptydrops.fdr = emptydrops.fdr, BPPARAM = cl, my.seed = my.seed)
  
  ### Computing basic metrics
  sobj <- QC.metrics(sobj = sobj, mt.genes.file = mt.genes.file, crb.genes.file = crb.genes.file, str.genes.file = str.genes.file, pcmito.range = pcmito.range, pcribo.range = pcribo.range, min.features = min.features, min.counts = min.counts, BPPARAM = cl)
}

### Filtering cells on metrics
sobj <- cells.QC.filter(sobj = sobj, pcmito.range = pcmito.range, pcribo.range = pcribo.range, min.features = min.features, min.counts = min.counts)

### Cell cycle prediction
sobj <- cell.cycle.predict(sobj = sobj, cc.pairs.file = cc.pairs.file, cc.genes.file = cc.genes.file, assay = 'RNA', nbin = 10, BPPARAM = cl)

### Tag control genes
sobj <- tag.ctrl.genes(sobj = sobj, ctrl.genes = ctrl.genes, ctrl.min.counts = 3, assay = "RNA")

### Filtering features (based on minimum cells covering)
sobj <- features.filter(sobj = sobj, min.cells = min.cells)

### Identification of doublets
sobj <- find.doublets(sobj = sobj, assay = assay, min.clust.size = min.clust.size)

### Building output directiry
out.dir <- paste(c(data.path, paste0('F', sobj@misc$params$min.features, '_C', sobj@misc$params$min.counts, '_M', paste(sobj@misc$params$pcmito.range, collapse = '-'), '_R', paste(sobj@misc$params$pcribo.range, collapse = '-'), collapse = '_'), 'DOUBLETSKEPT'), collapse = '/')
dir.create(path = out.dir, recursive = TRUE, showWarnings = TRUE)

### Saving non-normalized object
save(sobj, file = paste0(out.dir, '/', sample.name, '_DOUBLETSKEPT_NON-NORMALIZED.rda'), compress = "bzip2")

### Basic normalization
sobj <- sc.normalization(sobj = sobj, assay = assay, method = norm.method, features.n = 3000, vtr = NULL)

### Building normalized output dir
red.dir <- paste(out.dir, norm.method, sep = '/')
dir.create(path = red.dir, recursive = TRUE, showWarnings = TRUE)

### Saving normalized object
save(sobj, file = paste0(red.dir, '/', sample.name, '_DOUBLETSKEPT_', norm.method, '.rda'), compress = "bzip2")

if (tolower(norm.method) == 'sctransform') assay <- 'SCT'
### Dim reduction
sobj <- dimensions.reduction(sobj = sobj, reduction = reduction, assay = assay, max.dims = keep.dims)

### Building reduced output dir
clust.dir <- paste(red.dir, paste(c(reduction, keep.dims, keep.res), collapse = '_'), sep = '/')
dir.create(path = clust.dir, recursive = TRUE, showWarnings = TRUE)

## Clustering + uMAP
sobj <- louvain.cluster(sobj = sobj, reduction = reduction, assay = assay, max.dim = keep.dims, resolution = keep.res, out.dir = clust.dir, solo.pt.size = solo.pt.size)

## Technical plots
technical.plot(sobj = sobj, assay = assay, umap.name = paste(c(assay, reduction, keep.dims, 'umap'), collapse = '_'), ident = paste0(paste(c(assay, reduction, keep.dims), collapse = '.'), '_res.', keep.res), out.dir = clust.dir, multi.pt.size = multi.pt.size)

### Saving reduced object
save(sobj, file = paste0(clust.dir, '/', paste(c(sample.name, 'DOUBLETSKEPT', norm.method, reduction, vtr), collapse = '_'), '.rda'), compress = "bzip2")






## C. FILTERED ROLL, DOUBLETS REMOVED, WITH OT WITHOUT COVARIATE REGRESSION
##-------------------------------------------------------------------------

### PARAMS
assay <- 'RNA'
earlier.norm.method <- 'SCTransform' ## Allows to retrieve step B data but using another norm method for step C
norm.method <- 'SCTransform' ## or 'SCTransform'
reduction <- 'pca'
max.dims <- 80 ## depends on sample complexity and number of cells
vtr <- NULL


filtered.dir <- if (exists("sobj")) {
  paste(
    data.path,
    paste0(
      'F', sobj@misc$params$min.features,
      '_C', sobj@misc$params$min.counts,
      '_M', paste(sobj@misc$params$pcmito.range, collapse = '-'),
      '_R', paste(sobj@misc$params$pcribo.range, collapse = '-'),
      collapse = '_'
    ),
    sep = '/'
  )
} else {
  paste(
    data.path,
    paste0(
      'F', min.features,
      '_C', min.counts,
      '_M', paste(pcmito.range, collapse = '-'),
      '_R', paste(pcribo.range, collapse = '-'),
      collapse = '_'
    ),
    sep = '/'
  )
}

expected.rda <- paste0(
  paste(
    c(filtered.dir,
    'DOUBLETSKEPT',
    earlier.norm.method),
    collapse = '/'),
  '/',
  paste(c(sample.name, 'DOUBLETSKEPT', earlier.norm.method), collapse = '_'),
  '.rda'
)

### Creating parallel instance
cl <- create.parallel.instance(nthreads = nthreads)

if (file.exists(expected.rda)) {
  message(paste0("Loading ", expected.rda, " ..."))
  load(expected.rda)
} else {
  
  ### Loading raw count matrix
  sobj <- load.sc.data(data.path = data.path, sample.name = sample.name, droplets.limit = droplets.limit, emptydrops.fdr = emptydrops.fdr, BPPARAM = cl, my.seed = my.seed)
  
  ### Computing basic metrics
  sobj <- QC.metrics(sobj = sobj, mt.genes.file = mt.genes.file, crb.genes.file = crb.genes.file, str.genes.file = str.genes.file, pcmito.range = pcmito.range, pcribo.range = pcribo.range, min.features = min.features, min.counts = min.counts, BPPARAM = cl)
  
  ### Filtering cells on metrics
  sobj <- cells.QC.filter(sobj = sobj, pcmito.range = pcmito.range, pcribo.range = pcribo.range, min.features = min.features, min.counts = min.counts)
  
  ### Cell cycle prediction
  sobj <- cell.cycle.predict(sobj = sobj, cc.pairs.file = cc.pairs.file, cc.genes.file = cc.genes.file, assay = 'RNA', nbin = 10, BPPARAM = cl)
  
  ### Filtering features (based on minimum cells covering)
  sobj <- features.filter(sobj = sobj, min.cells = min.cells)
  
  ### Identification of doublets
  sobj <- find.doublets(sobj = sobj, assay = assay, min.clust.size = min.clust.size)
}

## Filter doublets (scds-hubrid + scDblFinder)
sobj <- filter.doublets(sobj = sobj, method = 'both')

## Tag control genes
sobj <- tag.ctrl.genes(sobj = sobj, ctrl.genes = ctrl.genes, ctrl.min.counts = 3, assay = "RNA")

### Building filtered output dir
filtered.dir <- paste(
  data.path,
  paste0(
    'F', sobj@misc$params$min.features,
    '_C', sobj@misc$params$min.counts,
    '_M', paste(sobj@misc$params$pcmito.range, collapse = '-'),
    '_R', paste(sobj@misc$params$pcribo.range, collapse = '-'),
    collapse = '_'
  ),
  sep = '/'
)
dir.create(filtered.dir, showWarnings = TRUE, recursive = TRUE)

### Saving normalized object
save(sobj, file = paste0(filtered.dir, '/', sample.name, '_FILTERED_NON-NORMALIZED.rda'), compress = "bzip2")

### Building normalized output dir
norm.dir <- paste(
  c(filtered.dir,
    if(is.null(vtr)) 'UNREGRESSED' else paste(vtr, collapse = '_'),
    norm.method),
  collapse = '/'
)
dir.create(path = norm.dir, recursive = TRUE, showWarnings = TRUE)

## Normalization + Dimension reduction
if (reduction %in% c("scbfa", "bpca")) {
  sobj <- binary.processing(sobj = sobj, assay = assay, method = reduction, norm = norm.method, features.n = 3000, max.dims = max.dims, vtr = vtr, vtr.scale = TRUE)
} else {
  sobj <- sc.normalization(sobj = sobj, assay = assay, method = norm.method, features.n = 3000, vtr = vtr)
  ### Saving normalized object
  save(sobj, file = paste0(norm.dir, '/', paste(c(sample.name, if(is.null(vtr)) 'UNREGRESSED' else paste(vtr, collapse = '_'), norm.method), collapse = '_'), '.rda'), compress = "bzip2")
  if (tolower(norm.method) == 'sctransform') assay <- 'SCT'
  sobj <- dimensions.reduction(sobj, reduction = reduction, assay = assay, max.dims = max.dims)
}

if (tolower(norm.method) == 'sctransform') assay <- 'SCT'

### Building reduced output dir
red.dir <- paste(
  c(norm.dir, reduction),
  collapse = '/'
)
dir.create(path = red.dir, recursive = TRUE, showWarnings = TRUE)

### Saving reduced object
save(sobj, file = paste0(red.dir, '/', paste(c(sample.name, if(is.null(vtr)) 'UNREGRESSED' else paste(vtr, collapse = '_'), norm.method, reduction), collapse = '_'), '.rda'), compress = "bzip2")

### QC histograms
QC.hist(sobj = sobj, out.dir = red.dir)

### Correlating reduction dimensions with biases and markers expression
dimensions.eval(sobj, reduction = reduction, assay = assay, eval.markers = eval.markers, out.dir = red.dir)

### Testing multiple clustering parameters (nb dims kept + Louvain resulution)
resvec <- seq(.1,1.2,.1)
dimsvec <- c(seq.int(5L, 49L, 2L),seq.int(50L,max.dims,5L))
clustering.eval.mt(sobj = sobj, reduction = reduction, assay = assay, dimsvec = dimsvec, resvec = resvec, out.dir = red.dir, solo.pt.size = solo.pt.size, BPPARAM = cl)

## PAUSE : Observe ALL UMAPS !

### Manually Choose optimal dims + resolution from clustree / uMAPS
keep.dims <- 27L
keep.res <- 0.8

### Building clustered output directory 
clust.dir <- paste(red.dir, paste(c(reduction, keep.dims, keep.res), collapse = '_'), sep = '/')
dir.create(path = clust.dir, recursive = TRUE, showWarnings = TRUE)

### Replotting final clusters
sobj <- louvain.cluster(sobj = sobj, reduction = reduction, assay = assay, max.dim = keep.dims, resolution = keep.res, out.dir = clust.dir, solo.pt.size = solo.pt.size)

### Technical plots
technical.plot(sobj = sobj, assay = assay, umap.name = paste(c(assay, reduction, keep.dims, 'umap'), collapse = '_'), ident = paste0(paste(c(assay, reduction, keep.dims), collapse = '.'), '_res.', keep.res), out.dir = clust.dir, multi.pt.size = multi.pt.size)

### Assessing clusters : finding markers
get.idents(sobj)
sobj <- find.markers.quick(sobj = sobj, assay = assay, ident = paste0(paste(c(assay, reduction, keep.dims), collapse = '.'), '_res.', keep.res), test.use = 'wilcox', min.pct = .75, logfc.threshold = .5, only.pos = TRUE, topn = 10, heatmap.cols = gradient.cols, out.dir = clust.dir)

### Assessing clusters : automatic cell type annotation
get.reductions(sobj)
sobj <- cells.annot(sobj = sobj, assay = assay, ident = paste0(paste(c(assay, reduction, keep.dims), collapse = '.'), '_res.', keep.res), umap.name = paste(c(assay, reduction, keep.dims, "umap"), collapse = "_"), singler.setnames = singler.setnames, clustifyr.setnames = clustifyr.setnames, cellassign.markers = markers.cellassign, sr.minscore = .25, cfr.minscore = .35, out.dir = clust.dir, solo.pt.size = solo.pt.size)

### Assessing clusters : Plotting control genes
sobj <- ctrl.umap.plot(sobj = sobj, ctrl.genes = ctrl.genes, assay = assay, umap.name = paste(c(assay, reduction, keep.dims, 'umap'), collapse = '_'), ident = paste0(paste(c(assay, reduction, keep.dims), collapse = '.'), '_res.', keep.res), out.dir = clust.dir, dimplot.cols = gradient.cols, multi.pt.size = 2)

### Assessing clusters : Plotting provided marker genes
sobj <- markers.umap.plot(sobj = sobj, markers = markers, assay = assay, umap.name = paste(c(assay, reduction, keep.dims, 'umap'), collapse = '_'), ident = paste0(paste(c(assay, reduction, keep.dims), collapse = '.'), '_res.', keep.res), out.dir = clust.dir, dimplot.cols = gradient.cols, multi.pt.size = 2)

### Saving final object
save(sobj, file = paste0(clust.dir, '/', paste(c(sample.name, if(is.null(vtr)) 'UNREGRESSED' else paste(vtr, collapse = '_'), norm.method, reduction, keep.dims, keep.res), collapse = '_'), '.rda'), compress = "bzip2")

### Building cerebro binary
sobj <- seurat2cerebro(sobj = sobj, assay = assay, ident = paste0(paste(c(assay, reduction, keep.dims), collapse = '.'), '_res.', keep.res), species = species.rdx, gmt.file = gmt.file, out.dir = clust.dir, nthreads = nthreads, only_pos = TRUE, min_pct = .75, thresh_logFC = .5, thresh_p_val = 1E-02, test = 'wilcox')

### Saving cerebro-compatible object
save(sobj, file = paste0(clust.dir, '/', paste(c(sample.name, if(is.null(vtr)) 'UNREGRESSED' else paste(vtr, collapse = '_'), norm.method, reduction, keep.dims, keep.res, 'CerebroReady'), collapse = '_'), '.rda'), compress = "bzip2")

## Count correction
source("SoupX.R")
# Load raw matrix from cell ranger
raw_counts <- Seurat::Read10X("/mnt/Data/Project/NF1-MPNST_Mouse/2020_02/counts/raw_feature_bc_matrix/")


path <- '/mnt/Data/Project/NF1-MPNST_Mouse/2020_02/'

# plot are stored in subfolder SOUPX in path
out <- remove_contamination(sobj, raw_counts,path, reduction = "SCT_pca_25_umap", reduction_keys = "sctpca25umap_", genes = "dtTomato")

### OTHERS : ADT

## ADT data (protein expression from CITE-seq)
ADT.assay <- 'ADT'
RNA.assay <- 'RNA'
norm.method <- 'LogNormalize'
gene.names <- c('CD3G', 'CD4', 'CTLA4', 'IL2RA', 'PDCD1', 'DPP4', 'MS4A1', 'CD24', 'SDC1', 'CR2', 'CD38', 'CD19')
GEdata <- '/home/job/WORKSPACE/SINGLECELL/RUN/IGR/191212_A00461_0092_AHH5H3DRXX/P30_FXDA/0732M_GE_cDNA/F200_C1000_M0-0.05_R0-1/UNREGRESSED_scBFA/scbfa_25_0.9/0732M_GE_cDNA_scbfa.rda'
cor.method <- 'spearman'
slot <- 'data'


cl <- create.parallel.instance(nthreads = nthreads)

### Loading raw count matrix
sobjADT <- load.sc.data(data.path = data.path, sample.name = sample.name, droplets.limit = 1E+6, emptydrops.fdr = NULL, BPPARAM = cl, my.seed = my.seed)

### Loading GE dataset
load(GEdata)

### Synching ADT to GE cells
adt.mtx <- sobjADT@assays[[RNA.assay]]@counts
rm(sobjADT)
adt.mtx <- adt.mtx[,colnames(adt.mtx) %in% colnames(sobj)]
sobj <- sobj[,colnames(sobj) %in% colnames(adt.mtx)]
if(!all(sort(colnames(sobj)) == sort(colnames(adt.mtx)))) message("SYNCH FAIL !")

### Merging
sobj[[ADT.assay]] <- Seurat::CreateAssayObject(adt.mtx)

### Normalization
sobj <- Seurat::NormalizeData(sobj, assay = assay, normalization.method = norm.method) ## Gave good VISUAL results. Avoid sctransfomr on such small dataset

### Computing correlations
cor.df <- data.frame(RNA_feature = gene.names, ADT_feature = rownames(sobj@assays[[ADT.assay]]@counts), stringsAsFactors = FALSE)

suppressWarnings(cor.unfiltered <- feature.cor(sobj = sobj, assay1 = RNA.assay, assay2 = ADT.assay, assay1.features = gene.names, assay2.features = rownames(sobj@assays[[ADT.assay]]@counts), slot = slot, cor.method = cor.method, zero.filter = FALSE))
suppressWarnings(cor.filtered <- feature.cor(sobj = sobj, assay1 = RNA.assay, assay2 = ADT.assay, assay1.features = gene.names, assay2.features = rownames(sobj@assays[[ADT.assay]]@counts), slot = slot, cor.method = cor.method, zero.filter = TRUE))
cor.df <- cbind(cor.df, cor.unfiltered, cor.filtered)

## Saving correlation results in @misc slot
sobj@assays[[ADT.assay]]@misc$cor <- cor.df

### Saving
sobj@assays[[RNA.assay]]@scale.data <- matrix()
sobj@assays[[ADT.assay]]@scale.data <- matrix()
save(sobj, file = paste0(data.path, '/', sample.name, '_GE+ADT_', norm.method, '.rda'), compress = "bzip2")




### OTHERS : TCR

## TCR enrichment data (quantification)
assay <- 'TCR'
norm.method <- 'LogNormalize'
reduction <- 'scbfa'
features.n <- 50 ## VERY important for TCR
max.dims <- 50

## Seurat object to insert TCR data. Will serve as basis for cells.
GEdata <- '/home/job/WORKSPACE/SINGLECELL/RUN/IGR/191212_A00461_0092_AHH5H3DRXX/P30_FXDA/0732M_ADT/0732M_ADT_GE+ADT_LogNormalize.rda'
file.exists(GEdata)


data.path <- dirname(GEdata)
cl <- create.parallel.instance(nthreads = nthreads)

### Loading raw count matrix
sobjTCR <- load.sc.data(data.path = data.path, sample.name = sample.name, assay = assay, droplets.limit = 1E+06, emptydrops.fdr = NULL, BPPARAM = cl, my.seed = my.seed)

### Loading GE dataset
load(GEdata)

### Synching Seurat objects
## We need to filter cells to community
sobjTCR <- sobjTCR[rownames(sobjTCR) %in% rownames(sobj), colnames(sobjTCR) %in% colnames(sobj)]
sobj <- sobj[rownames(sobj) %in% rownames(sobjTCR), colnames(sobj) %in% colnames(sobjTCR)]
if(!all(rownames(sobj) == rownames(sobjTCR))) message("SYNCH FAIL on rows !")
if(!all(colnames(sobj) == colnames(sobjTCR))) message("SYNCH FAIL on columns !")

### Merging
sobj[[assay]] <- Seurat::CreateAssayObject(sobjTCR@assays[[assay]]@counts)
rm(sobjTCR)

### Building reduction data dir
red.dir <- paste(c(data.path, features.n, norm.method, reduction), collapse = '/')
dir.create(red.dir, recursive = TRUE)

### Normalization / reduction
if (reduction %in% c('scbfa', 'bpca')) {
  sobj <- binary.processing(sobj = sobj, assay = 'TCR', norm = norm.method, method = reduction, features.n = features.n, max.dims = max.dims, vtr = NULL)
} else {
  sobj <- sc.normalization(sobj = sobj, assay = assay, method = norm.method, features.n = features.n, vtr = NULL)
  sobj <- Seurat::ScaleData(sobj, assay = assay, features = rownames(sobj@assays[[assay]]@counts))
  sobj <- dimensions.reduction(sobj, reduction = reduction, assay = assay, max.dims = max.dims)
}

### Evaluating clustering parameters
resvec <- seq(.1,1.5,.1)
dimsvec <- seq.int(5L, max.dims, 2L)
clustering.eval.mt(sobj = sobj, reduction = reduction, assay = assay, dimsvec = dimsvec, resvec = resvec, out.dir = red.dir, solo.pt.size = solo.pt.size, BPPARAM = cl)

### Building clustering output dir

### Replotting final clusters
keep.dims <- 20L
keep.res <- .8
clust.dir <- paste0(red.dir, '/', paste(c(reduction, keep.dims, keep.res), collapse = '_'))
dir.create(clust.dir, recursive = TRUE)
sobj <- louvain.cluster(sobj = sobj, reduction = reduction, assay = assay, max.dim = keep.dims, resolution = keep.res, out.dir = clust.dir, solo.pt.size = solo.pt.size)

### Rough markers identification
sobj <- find.markers.quick(sobj = sobj, assay = assay, ident = paste0(assay, '.', reduction, '.', keep.dims, '_res.', keep.res), test.use = 'wilcox', min.pct = .75, logfc.threshold = .5, only.pos = TRUE, topn = 10, heatmap.cols = gradient.cols, out.dir = clust.dir)

### Saving
save(sobj, file = paste0(clust.dir, '/', sample.name, '_GE+ADT+TCR_', norm.method, '.rda'), compress = "bzip2")

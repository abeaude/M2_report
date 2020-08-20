require(furrr)
options(future.globals.maxSize = 10^10)
plan(multicore,workers = 6)

object_name <- "sobj"
assay <- "SCT"
reference <- NULL
file_list <- c('/mnt/Data/Project/NF1-MPNST/NF1-MPNST_20181205/ANALYSIS/DEFAULT/NF1-MPNST_20181205_SCT_DEFAULT_pca70_resolution0.7.rda','/mnt/Data/Project/NF1-MPNST/NF1-MPNST_20190227/ANALYSIS/DEFAULT/NF1-MPNST_20190227_SCT_DEFAULT_pca50_resolution0.6.rda','/mnt/Data/Project/NF1-MPNST/NF1-MPNST_20190731/ANALYSIS/DEFAULT/NF1-MPNST_20190731_SCT_DEFAULT_pca40.rda')

assay_vec <- rep(assay,length(sobj_list))


integration.metrics <- function(sobj.list,sobj.features,k.anchor,k.filter,k.score,max.features,nn.method,dims,assay_vec,reference,n){
  print(n)
  sobj.anchors <- Seurat::FindIntegrationAnchors(object.list = sobj_list, assay = assay_vec,reference = reference, normalization.method = "SCT", anchor.features = sobj.features, verbose = FALSE, scale = FALSE, reduction = 'cca',k.anchor = k.anchor,k.filter = k.filter,k.score = k.score,max.features = max.features,nn.method = nn.method,dims = 1:dims,l2.norm = TRUE)
  
  suppressWarnings(sobj_integrated <- Seurat::IntegrateData(anchorset = sobj.anchors, normalization.method = "SCT", verbose = FALSE))
  sobj_integrated <- dimensions.reduction(sobj = sobj_integrated, reduction = 'pca', assay = 'integrated', vtr = NULL, red.dims = 1:100, my.seed = 1337)
  rm(sobj.anchors)
  gc(verbose = FALSE)
  max.k <- 300 
  ls <- Seurat::LocalStruct(sobj_integrated,grouping.var = "orig.ident",reduction = 'pca',reduced.dims = 1:dims, orig.dims = 1:dims, verbose = FALSE)
  mm <- max.k - Seurat::MixingMetric(sobj_integrated,grouping.var = "orig.ident",reduction = 'pca',dims = 1:dims,max.k = max.k, k = k.anchor, verbose = FALSE)
  rm(sobj_integrated)
  gc(verbose = FALSE)
  return(list(mean(unname(obj = unlist(x = ls)),na.rm = TRUE),mean(mm,na.rm = TRUE)))
}

k.anchor <- seq.int(5,10,1)
k.filter <- 200#seq.int(100,500,50)
k.score <- 30# seq.int(10,100,10)
max.features <- 200#seq.int(100,500,100)
nn.method <- c("rann","annoy")
dims <- c(seq.int(10,50,10))
library(magrittr)
options(future.globals.maxSize = 4000 * 1024^2)
message("Processing step")
start_time <- Sys.time()
params <- tidyr::expand_grid(nn.method,dims,max.features,k.score,k.filter,k.anchor) %>% dplyr::mutate(n = seq.int(1,dim(.)[1]))
print(dim(params)[1])
results <-  params %>% furrr::future_pmap(.f=integration.metrics, sobj.features = sobj.features, sobj.list = sobj_list, assay_vec = assay_vec,reference = reference) %>% dplyr::mutate(local_structure = purrr::map_dbl(results,1), mixing_metric = purrr::map_dbl(results,2)) %>% dplyr::select(-results)
end_time <- Sys.time()
message('DONE')

save(start_time,end_time, params,results,file = "test_benchmarck_integration.rda")




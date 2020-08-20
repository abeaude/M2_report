liger_integration <- function(sobj_list, k = 20, lambda = 5, split.by = "orig.ident", project = 'LIGER_integration', verbose = FALSE, seed = NULL){
  assertthat::assert_that(are_seurat(sobj_list))
  assertthat::assert_that(is.numeric(k), is.numeric(lambda))
  assertthat::assert_that(is.character(project))
  
  seed <- seed %||% sobj_list[[1]]@misc$params$seed %||% 1337L
  
  sobj <- merge(sobj_list[[1]], sobj_list[2:length(sobj_list)], merge.data = FALSE, project = project)
  assertthat::assert_that(has_metadata(sobj, split.by))
  
  Seurat::DefaultAssay(sobj) <- "RNA"

  sobj <- sobj %>%
    Seurat::NormalizeData(verbose = verbose) %>%
    Seurat::FindVariableFeatures(verbose = verbose)  %>%
    Seurat::ScaleData(split.by = split.by, do.center = FALSE, verbose = verbose)  %>%
    SeuratWrappers::RunOptimizeALS(k = k, lambda = lambda, split.by = split.by, rand.seed = seed) %>%
    SeuratWrappers::RunQuantileNorm(split.by = split.by, reduction.name = "RNA_iNMF", reduction.key = "RNAiNMF_")
  
  message("To perform any downstream analysis you can use the reduced space with the key : 'iNMF'\n",
          ncol(Seurat::Reductions(sobj, "RNA_iNMF")), " available dimensions")
  sobj
}

create_liger_object <- function(sobj_list){
  assertthat::assert_that(are_seurat(sobj_list), is_named(sobj_list))
  sobj_list %>%
    purrr::map(~ .@assays$RNA@counts) %>% 
    liger::createLiger()
}



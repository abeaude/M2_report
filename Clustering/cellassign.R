cellassign.by.cluster <- function(sobj, metadata.cellassign = 'cellassign'){
  if(is.null(sobj)){
    stop("No seurat object provided")
  }
  cluster_assign <- tibble::tibble(cluster = as.numeric(unique(sobj[['seurat_clusters']])[[1]])-1,
                                   assign = "Unknown")
  table_cell <- tibble::tibble(celltype=sobj[[metadata.cellassign]],barcode = rownames(sobj[[]]), cluster=as.numeric(sobj[['seurat_clusters']][[1]])-1)
  
  for(clust in cluster_assign$cluster) {
    freqs <- table_cell %>% dplyr::filter(cluster == clust) %$%
      table(celltype) %>% tibble::as_tibble() 
    
    num_unk <- dplyr::filter(freqs,celltype =='other')
    freqs <- dplyr::filter(freqs,celltype !='other')
    putative_type <- freqs %>% 
      dplyr::top_n(n = 1,n) %>% 
      dplyr::select(celltype) %>% 
      dplyr::pull()
    cluster_assign %<>% dplyr::mutate(assign = ifelse(cluster == clust,putative_type,assign))
  }
  sobj[["cellassign_cluster"]] <- dplyr::full_join(table_cell,cluster_assign,by = 'cluster') %$% assign
  return(sobj)
}

cell.annotation.cellassign <- function(sobj, markers, GPU = FALSE, by.cluster = TRUE, seed = NULL, verbose = FALSE, min_delta = 2, learning_rate = 1e-2, shrinkage = TRUE, num_runs = 1){
  assertthat::assert_that(is_seurat(sobj))
  assertthat::assert_that(is.character(markers), is_named(markers))
  seed <- seed %||% sobj@misc$params$seed %||% 1337L
  if(!GPU){
    message("Tensorflow will use the CPU")
    Sys.setenv("CUDA_DEVICE_ORDER"="PCI_BUS_ID")
    Sys.setenv("CUDA_VISIBLE_DEVICES" = "-1")
  }
  
  not.in.sobj <- which(!markers %in% rownames(sobj@assays$RNA@counts))
  if(any(not.in.sobj)){
    markers <- markers[-not.in.sobj]
  }
  marker_mat <- cellassign::marker_list_to_mat(split(unname(markers),names(markers)))
  #From cellassign help
  # The assay from the input#' SingleCellExperiment to use: this assay should always represent raw counts.
  s <- Seurat::as.SingleCellExperiment(sobj,assay = "RNA") %>%
    scran::computeSumFactors() %>%
    SingleCellExperiment::sizeFactors()
  
  # set.seed to reproductibility
  set.seed(seed)
  fit <- cellassign::cellassign(exprs_obj = Seurat::as.SingleCellExperiment(sobj,assay = "RNA")[rownames(marker_mat),], 
                                marker_gene_info = marker_mat, 
                                s = s, 
                                min_delta = min_delta,
                                learning_rate = learning_rate, 
                                shrinkage = shrinkage,
                                verbose = verbose,
                                num_runs = num_runs)
  sobj[["cellassign"]] <- fit$cell_type
  sobj@misc$cellassign_lls <- fit$lls
  sobj@misc$params$cellassign <- fit$mle_params
  if(by.cluster){
    sobj <- cellassign.by.cluster(sobj = sobj)
  }
  return(sobj)
}

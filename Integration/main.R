integration_metrics <- function(sobj_integrated, group_by = "orig.ident", reduction = 'pca', max.k, k.anchor, dims){
  ls <- Seurat::LocalStruct(sobj_integrated,grouping.var = group_by, reduction = reduction, reduced.dims = 1:dims, orig.dims = 1:dims, verbose = FALSE) %>%
    unlist() %>%
    unname() %>%
    mean(na.rm = TRUE)
  
  mm <- max.k - Seurat::MixingMetric(sobj_integrated,grouping.var = group_by, reduction = reduction, dims = 1:dims, max.k = max.k, k = k.anchor, verbose = FALSE) %>%
    mean(na.rm = TRUE)
      
  return(list("Local Strucutre" = ls,
              "Mixing Metric" = mm))
}

plot_integration_metrics <- function(metrics_integration){
  
}

create_sobj_list <- function(...){
  paths <- rlang::list2(...) %>% unlist()
  assertthat::assert_that(is.character(paths))
  purrr::walk(paths, ~assertthat::assert_that(assertthat::is.readable(.)))
  
  sobj_list <- list()
  for(i in seq_along(paths)){
    if(fs::path_ext(paths[i]) == 'rda' || fs::path_ext(paths[i]) == 'RData'){
      sobj_list[[i]] <- get(load(paths[i]))
    } else if (fs::path_ext(path[i]) == 'rds'){
      sobj_list[[i]] <- readRDS(paths[i])
    } else {
      warning("Unknown extension : '",fs::path_ext(paths[i]), "'. Not added to the list")
    }
  }
  # Add names according to project name of each object
  purrr::set_names(sobj_list, purrr::map(sobj_list, ~ .x@project.name)) %>%
  # MAke cells unique
  unique_cell_names()
}

unique_cell_names <- function(sobj_list){
  initial_order <- names(sobj_list)
  unique_names <- purrr::map(sobj_list, colnames) %>% 
    unlist(use.names = FALSE) %>% 
    purrr::set_names(rep(names(sobj_list), times = purrr::map(sobj_list, ~ colnames(.) %>% length()))) %>%
    vctrs::vec_as_names(repair = "unique") %>%
    {split(unname(.), names(.))}
  
  unique_names = unique_names[order(names(unique_names))]
  sobj_list = sobj_list[order(names(sobj_list))]
  
  purrr::map2(sobj_list, unique_names, ~ Seurat::RenameCells(.x,new.names = .y)) %>%
    magrittr::extract(initial_order)
}

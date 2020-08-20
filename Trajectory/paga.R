ti_paga <- function(params,dir_name, path, filename){
  # https://rdrr.io/github/milescsmith/ReductionWrappers/src/R/paga.R
  # https://romanhaa.github.io/blog/paga_to_r/
  
  # adata <- sc$read_loom(filename = fs::path(path,paste0(filename,".loom")), sparse = TRUE)
  adata <- sc$read(filename = fs::path(path,paste0(filename,".h5ad")))
  if(params[['use_rep']] == "X_pca") sc$tl$pca(adata,n_comps = params[['n_comps']], svd_solver = 'arpack')
  if(params[['connectivities']] == "scanpy"){
    sc$pp$neighbors(adata, random_state = params[['random_state']], use_rep = params[['use_rep']])
  } else if(params[['connectivities']] == "seurat"){
    connectivities <- readRDS(fs::path(path, paste0(filename,'_connectivities.rds')))
    
    adata$obsp <- reticulate::dict("connectivities" = connectivities)
    adata$uns$neighbors <- reticulate::dict('connectivities_key'= 'connectivities')
    # if params is absent will compute n_neighbors automatically 
    
    # # paga V1.0
    # adata$obsp['connectivities'] <- connectivities
    # # paga v1.2
    # adata$obsp['distances'] <- distances
    
    # For the moment with reticulate we can only access the global environment inside a python session
    # and setting obsp from the r session doesn't work. 
    # Copy from local env to global ? 
    
    # reticulate::py_run_string("r.adata.obsp['connectivities'] = r.connectivities")
  }
  # denoise the graph by recomputing it in the first few diffusion components
  if (params[['denoise_graph']]) {
    sc$tl$diffmap(adata, n_comps=params[['n_comps_diff']])
    sc$pp$neighbors(adata, random_state = params[['random_state']], use_rep = 'X_diffmap')
  }
  
  if(params[['groups']] == "leiden") {
    sc$tl$leiden(adata, resolution = params[['resolution']], random_state = params[['random_state']])
    adata$obs[['leiden']] <- adata$obs[['leiden']]$values$astype('int') %>% 
      {reticulate::py_to_r(.) + 1L} %>% 
      as.integer() %>%
      reticulate::r_to_py()
    
    adata$obs[['leiden']] <- adata$obs[['leiden']]$astype('category')
  }
    
  
  sc$tl$paga(adata, groups = params[['groups']], model = params[['paga_version']])

  if (!params[['denoise_graph']]) sc$tl$diffmap(adata, n_comps=params[['n_comps_diff']])
  
  # Pseudotime need a root cell 
  # paste0("r.adata.uns['iroot'] = int(r.np.where(r.adata.obs.index == '",params$start_cell,"')[0][0])") %>%
  #   reticulate::py_run_string()
  # python indices start at 0 not 1
  # adata$uns$iroot <- which(adata$obs$index$values %>% reticulate::py_to_r() == params[['start_cell']]) %>% as.integer() - 1L
  sc$tl$dpt(adata, n_dcs = params[['n_dcs']]) # min(adata$obsm['X_diffmap']$shape[1] %>% reticulate::py_to_r() , 10L)
  
  
  plot_paga <- sc$pl$paga(adata, show = FALSE, threshold = params[["connectivities_cutoff"]], edge_width_scale = 0.5, node_size_scale = 2)
  plt$savefig(fs::path(path,dir_name,paste0(filename,"_paga_plot.png")), dpi = 150L, bbox_inches = 'tight')
  sc$tl$draw_graph(adata, random_state = params[['random_state']], init_pos = 'paga')
  
  #saving connectivities cutoff in adata object
  adata$uns$update("connectivities_cutoff" = params[["connectivities_cutoff"]])
  # adata$uns['user_params'] <- reticulate::dict("connectivities_cutoff" = params[["connectivities_cutoff"]])
  
 return(adata)
}

run_paga <- function(path,filename){
  path %<>% fs::path("paga")
  paga_params <- yaml::read_yaml(fs::path(path,"paga_params.yaml")) %>%
    filter_computed_traj(path)
  
  create_trajectory_dir(paga_params,path)
  
  # conda_env <- "scanpy"
  # # Test if environment exist with the required packages  
  # if(!any(reticulate::conda_list() %>% dplyr::pull(name) == conda_env)){
  #   message("Installing...")
  #   reticulate::conda_install(envname = conda_env, packages = c("leidenalg","scanpy", "loompy"), channel = c("conda-forge", "bioconda"))
  #   reticulate::conda_install(envname = conda_env, packages = c("fa2"), pip = TRUE)
  # }
  # reticulate::use_condaenv(condaenv = conda_env)  
  sc <<- reticulate::import("scanpy", convert = FALSE)
  pd <<- reticulate::import("pandas", convert = FALSE)
  np <<- reticulate::import("numpy", convert = FALSE)
  plt <<- reticulate::import("matplotlib.pyplot", convert = FALSE)
  on.exit(rm(np,plt,sc,pd, envir = rlang::global_env()))
  
  output <- purrr::imap(paga_params,ti_paga, path = path, filename = filename) %>%
    purrr::map(wrapper_paga)
  
  
  return(output)
}

wrapper_paga <- function(adata){
  dimred_name <- paste0("X_draw_graph_",adata$uns['draw_graph']['params']['layout'] %>% reticulate::py_to_r())
  connectivity_cutoff <- adata$uns['connectivities_cutoff'] %>% reticulate::py_to_r()
  
  # reticulate::py_run_string("r.grouping = r.adata.obs[r.adata.uns['paga']['groups']].astype(str)")
  group <- adata$uns['paga']['groups'] %>% reticulate::py_to_r()
  grouping <- adata$obs[group]$astype('str') %>% 
    reticulate::py_to_r() %>%
    tibble::enframe(name = "cell_id", value = "group_id")
  # grouping <- adata$obs[[adata$uns[['paga']][['groups']]]] %>% 
  #   reticulate::py_to_r() %>% 
  #   tibble::enframe(name = "cell_id", value = "group_id")
  root <- adata$obs$index[adata$uns[['iroot']]] %>% reticulate::py_to_r()
  
  # reticulate::py_run_string('r.millestone = r.pd.DataFrame(r.np.triu(r.adata.uns["paga"]["connectivities"].todense(), k = 0), 
  #                                   index = r.adata.obs.leiden.cat.categories, 
  #                                   columns = r.adata.obs.leiden.cat.categories).stack().reset_index()')
  
  milestone_network <- adata$uns["paga"]["connectivities"]$todense() %>% 
    np$triu(k=0) %>%
    pd$DataFrame(index = adata$obs[group]$cat$categories, columns = adata$obs[group]$cat$categories) %>%
    reticulate::py_to_r() %>%
    tibble::as_tibble(rownames = "from") %>%
    tidyr::pivot_longer(cols = -from, names_to = 'to', values_to = "length") %>%
    dplyr::filter(length >= connectivity_cutoff) %>%
    dplyr::mutate(directed = FALSE, to = stringr::str_trim(to,side = "right"))
  
  millestones <- dplyr::select(milestone_network, from, to) %>% unlist(use.names = FALSE) %>% unique()
  # millestones <- unique(grouping[['group_id']])
  # reticulate::py_run_string(paste0("r.dim_red = r.pd.DataFrame([x for x in r.adata.obsm['",dimred_name, "'].T]).T"))
  # dimred <- dim_red %>% magrittr::set_colnames(stringr::str_c("comp_",1:ncol(.))) %>%
  #   tibble::add_column(cell_id = grouping$cell_id)
  
  dimred <- # check with higher dim ?
    adata$obsm[dimred_name] %>% 
    reticulate::py_to_r()  %>%
    magrittr::set_colnames(stringr::str_c("comp_",1:ncol(.))) %>%
    tibble::as_tibble() %>%
    tibble::add_column(cell_id = grouping$cell_id)
  
  dimred_milestones <- dimred %>%
    dplyr::inner_join(grouping, by = 'cell_id') %>%
    dplyr::group_by(group_id) %>%
    dplyr::summarise(comp_1 = mean(comp_1), comp_2 = mean(comp_2), .groups = "drop") %>%
    dplyr::filter(group_id %in% millestones) %>%
    tibble::column_to_rownames("group_id") %>%
    as.matrix()
  
  dimred <- dimred %>%
    tibble::column_to_rownames("cell_id") %>%
    as.matrix()
  
  # dimred_milestones <- adata$uns$paga$pos %>% 
  #   reticulate::py_to_r() %>% 
  #   magrittr::set_colnames(stringr::str_c("comp_",1:ncol(.))) %>% 
  #   magrittr::set_rownames(unique(grouping$group_id) %>% stringr::str_sort(numeric = TRUE))
  
  # branch_progressions <- adata$obs %>%
  #   reticulate::py_to_r() %>%
  #   tibble::rownames_to_column("cell_id") %>%
  #   dplyr::select(dpt_pseudotime, leiden, cell_id) %>%
  #   dplyr::mutate(dpt_pseudotime = ifelse(is.infinite(dpt_pseudotime), 1, dpt_pseudotime), branch_id = as.character(leiden)) %>%
  #   dplyr::group_by(leiden) %>%
  #   dplyr::mutate(percentage = (dpt_pseudotime - min(dpt_pseudotime)) / (max(dpt_pseudotime) - min(dpt_pseudotime))) %>%
  #   dplyr::ungroup(leiden) %>%
  #   dplyr::select(cell_id,branch_id,percentage)
  # 
  # branches <- adata$obs %>%
  #   reticulate::py_to_r() %>%
  #   dplyr::select(leiden,dpt_pseudotime) %>%
  #   dplyr::group_by(leiden) %>%
  #   dplyr::mutate(length = max(dpt_pseudotime) - min(dpt_pseudotime), branch_id = leiden) %>%
  #   dplyr::ungroup(leiden) %>%
  #   dplyr::select(branch_id, length) %>%
  #   dplyr::distinct() %>%
  #   dplyr::arrange(branch_id) %>%
  #   dplyr::mutate(branch_id = as.character(branch_id), directed = TRUE)
  # 
  # branch_network <- dplyr::select(milestone_network,from,to)
  # average_pseudotime <- adata$obs %>% 
  #   reticulate::py_to_r() %>% 
  #   dplyr::select(leiden,dpt_pseudotime) %>%
  #   dplyr::group_by(leiden) %>%
  #   dplyr::summarise(mean = mean(dpt_pseudotime))
  # 
  # branch_network <- dplyr::inner_join(branch_network, average_pseudotime, by = c("from" = "leiden")) %>%
  #   dplyr::inner_join(average_pseudotime, by = c('to' = "leiden")) %>%
  #   dplyr::rename(from_dpt = mean.x, to_dpt = mean.y) %>%
  #   dplyr::mutate(from_1 = ifelse(from_dpt > to_dpt, to, from), to_1 = ifelse(from_dpt > to_dpt, from, to), from = from_1, to = to_1) %>%
  #   dplyr::select(from,to)
  
  # grouping : adata$uns[['paga']][['groups']]
  # root_cell : adata$obs$index[adata$uns[['iroot']]]
  
  expression <- adata$layers['norm_data'] %>% 
    reticulate::py_to_r() %>%
    magrittr::set_colnames(adata$var_names$values %>% reticulate::py_to_r() %>% as.character()) %>%
    magrittr::set_rownames(adata$obs_names$values %>% reticulate::py_to_r() %>% as.character()) %>%
    as.matrix() %>% 
    as("dgCMatrix")
  
  count <- adata$X %>% 
    reticulate::py_to_r() %>%
    magrittr::set_colnames(adata$var_names$values %>% reticulate::py_to_r() %>% as.character()) %>%
    magrittr::set_rownames(adata$obs_names$values %>% reticulate::py_to_r() %>% as.character()) %>%
    as.matrix() %>% 
    as("dgCMatrix")
  
  out <- dynwrap::wrap_data(cell_ids = adata$obs_names$values %>%  reticulate::py_to_r() %>% as.character()) %>% 
    dynwrap::add_dimred_projection(
      grouping = grouping, 
      milestone_network = milestone_network,
      dimred = dimred,
      dimred_milestones = dimred_milestones
    ) %>% dynwrap::add_root(root_cell_id = root) %>%
    dynwrap::add_expression(counts = count,
                            expression = expression) %>%
    dynwrap::add_pseudotime()
  
  return(out)
}

data_prep_paga <- function(sobj, filename, path, reductions, cell_identities, seed, start_cell){
  fs::path(path,"paga") %>% fs::dir_create()
  # https://github.com/mojaveazure/loomR/issues/40
  # If there are any NA in acharacter columns conversion will not work
  # NA are mainly in columns from cell annotation 
  # easier to remove it for now 
  # cols_to_remove <- c("SR_MouseRNAseqData_cells","SR_MouseRNAseqData_clust","SR_ImmGenData_cells","SR_ImmGenData_clust","CFR_ref_MCA_cells","CFR_ref_MCA_clust","CFR_ref_tabula_muris_drop_cells","CFR_ref_tabula_muris_drop_clust","CFR_ref_tabula_muris_facs_cells","CFR_ref_tabula_muris_facs_clust","CFR_ref_moca_main_cells","CFR_ref_moca_main_clust","CFR_ref_immgen_cells","CFR_ref_immgen_clust","CFR_ref_mouse.rnaseq_cells","CFR_ref_mouse.rnaseq_clust")
  paga_params_default <- list("n_comps" = 50L,
                              "n_neighbors" = 15L,
                              "denoise_graph" = FALSE,
                              "start_cell" = start_cell,
                              "use_rep" = 'X_pca',
                              "random_state" = seed,
                              "resolution" = 1,
                              "groups" = "leiden",
                              "n_dcs" = 10L,
                              "n_comps_diff" = 15L,
                              "connectivities" = "scanpy",
                              "paga_version" = 'v1.2', 
                              "connectivities_cutoff" = 0.01)

  cols_to_remove <- Filter(function(x) x,colnames(sobj@meta.data) %>% purrr::set_names(.) %>% purrr::map_lgl( ~ any(is.na(sobj[[.]])))) %>% names()
  if(length(cols_to_remove) > 0){
    warning("(PAGA) \nDuring the conversion in loom file the following columns : \n" ,stringr::str_c("\t- ",cols_to_remove, collapse = "\n"), " \nwere removed before converting to loom because of the presence of NA", call. = FALSE)
  }

  # https://github.com/satijalab/seurat/issues/2017
  # sobj[[]] <- NULL is only upoorted for one columns at the moment
  for(i in cols_to_remove){
    sobj[[i]] <- NULL
  }

  if(any(cols_to_remove == cell_identities)){
    warning(cell_identities, " was removed from the metadata object because of the presence of NA, switching to seurat_cluster instead")
    cell_identities <- 'seurat_cluster'
  }
  ## Correct encoding error to be able to import the file in scanpy
  colnames(sobj[[]]) %>% 
    purrr::set_names(.,.) %>% 
    purrr::map(~ sobj@meta.data[[.]] %>% tools::showNonASCII()) %>% 
    purrr::compact() %>%
    purrr::imap(function(.y,col_name){
      x <- stringi::stri_trans_general(sobj@meta.data[[col_name]],"latin-ascii") 
      x <- iconv(as.character(x), "", "ASCII", "byte")
      Encoding(x) <- "latin1"
      x
    }) %>%
    purrr::iwalk(function(data,col){
      sobj@meta.data[col] <<- data
    })
  
  # available.reductions <- purrr::map_chr(reductions,~stringr::str_extract(names(sobj@reductions), .x) %>% purrr::discard(is.na))
  # Only PCA is useful with PAGA
  # filter reductions to extract PCA (not umap in the name)
  pca_red_name <- reductions %>% 
    stringr::str_subset("(?i)umap", negate = TRUE) %>%
    stringr::str_subset("(?i)pca") %>%
    magrittr::extract(1) %>% 
    stringr::str_remove_all("\\(\\?i\\)\\^|\\$") %>%
    stringr::str_c("cell", "embeddings", sep = "_")
  
  paga_params <- list('default' = paga_params_default)
  
  paga_params_default[['use_rep']] <- pca_red_name
  paga_params[['seurat_pca']] <- paga_params_default
  
  paga_params_default[['groups']] <- cell_identities
  paga_params[['seurat_pca_clustering']] <- paga_params_default
  
  # only if connectivities is available in the object
  # if(length(sobj@graphs) != 0){
  #   paga_params_default[['connectivities']] <- "seurat"
  #   paga_params_default[['paga_version']] <- "v1.0"
  #   paga_params[["seurat_pca_clustering_connectivities"]] <- paga_params_default
  #   saveRDS(sobj@graphs$SCT_snn, file = fs::path(path,'paga', paste0(filename,'_connectivities.rds')))
  # }
  
  
  yaml::write_yaml(paga_params, file = fs::path(path,"paga","paga_params.yaml"))
  pa_sobjSCT.loom <- suppressMessages(Seurat::as.loom(sobj, filename = fs::path(path,'paga', paste0(filename,'.loom')), verbose = FALSE))
  # remove loom file
  pa_sobjSCT.loom$close_all()
  # once it's converted to loom file convert it to h5ad from scanpy (save space and time)
  # conda_env <- "scanpy"
  # #Test if environment exist with the required packages
  # if(!any(reticulate::conda_list() %>% dplyr::pull(name) == conda_env)){
  #   message("Installing...")
  #   reticulate::conda_install(envname = conda_env, packages = c("leidenalg","scanpy"), channel = c("conda-forge", "bioconda"))
  #   reticulate::conda_install(envname = conda_env, packages = c("fa2", "loompy"), pip = TRUE)
  # }
  # reticulate::use_condaenv(condaenv = conda_env)
  sc <- reticulate::import("scanpy", convert = FALSE)
  adata <- sc$read_loom(filename = fs::path(path, "paga",paste0(filename,".loom")), sparse = TRUE)
  adata$uns <- reticulate::dict(iroot = which(adata$obs$index$values %>% reticulate::py_to_r() == start_cell) %>% as.integer() - 1L)
  adata$write(filename = fs::path(path, "paga",paste0(filename,".h5ad")))

  fs::file_delete(fs::path(path, "paga",paste0(filename,".loom")))
  
  # adata <- seurat_to_anndata(sobj)
  # adata$write(fs::path(path, "paga",paste0(filename,".h5ad")))
  # write.table(connectivities, fs::path(path,'paga', paste0(filename,'_connectivities.csv')), row.names = FALSE, col.names = FALSE, sep = ";", quote = FALSE)
}

# # How to update the milestone network by filtering some edges based on connectivity cutoff
# dynplot::plot_dimred(trajectory, color_cells = "grouping", label_milestones = TRUE)
# trajectory$milestone_network %<>%
#   dplyr::filter(length > 0.1)
# dynplot::plot_dimred(trajectory, color_cells = "grouping", label_milestones = TRUE)

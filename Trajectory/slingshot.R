wrapper_sligshot <- function(sds,dimred_name, groups){
  # Adapted from https://github.com/dynverse/ti_slingshot/blob/master/package/R/ti_slingshot.R
  
  #juste après le slingshot:
  # start_cell <- apply(slingshot::slingPseudotime(sds), 1, min) %>% sort() %>% head(1) %>% names()
  # start.clus <- clusters[[start_cell]]
  
  #   ____________________________________________________________________________
  #   Create output                                                           ####
  
  # collect milestone network
  lineages <- slingshot::slingLineages(sds)
  lineage_ctrl <- slingshot::slingParams(sds)
  
  cluster_network <- lineages %>%
    purrr::map_df(~ tibble::tibble(from = .[-length(.)], to = .[-1])) %>%
    unique() %>%
    dplyr::mutate(
      length = lineage_ctrl$dist[cbind(from, to)],
      directed = TRUE
    )
  
  # collect dimred
  # dimred <- SingleCellExperiment::reducedDim(sds)
  # collect clusters
  cluster <- slingshot::slingClusterLabels(sds)
  
  # collect progressions
  adj <- slingshot::slingAdjacency(sds)
  lin_assign <- apply(slingshot::slingCurveWeights(sds), 1, which.max)
  
  progressions <- purrr::map_df(seq_along(lineages), function(l) {
    ind <- lin_assign == l
    lin <- lineages[[l]]
    pst.full <- slingshot::slingPseudotime(sds@int_metadata$slingshot, na = FALSE)[,l]
    pst <- pst.full[ind]
    means <- sapply(lin, function(clID){
      stats::weighted.mean(pst.full, cluster[,clID])
    })
    non_ends <- means[-c(1,length(means))]
    edgeID.l <- as.numeric(cut(pst, breaks = c(-Inf, non_ends, Inf)))
    from.l <- lineages[[l]][edgeID.l]
    to.l <- lineages[[l]][edgeID.l + 1]
    m.from <- means[from.l]
    m.to <- means[to.l]
    
    pct <- (pst - m.from) / (m.to - m.from)
    pct[pct < 0] <- 0
    pct[pct > 1] <- 1
    
    a <- tibble::tibble(cell_id = names(which(ind)), from = from.l, to = to.l, percentage = pct)
  })
  
  grouping <- SingleCellExperiment::colData(sds)[groups] %>% 
    tibble::as_tibble(,rownames = 'cell_id') %>% 
    dplyr::rename(group_id = {{groups}})
  #   ____________________________________________________________________________
  #   Save output                                                             ####
  output <-  dynwrap::wrap_data(cell_ids = colnames(sds)) %>% 
    dynwrap::add_expression(counts = Matrix::t(SingleCellExperiment::counts(sds)) ,
                            expression = Matrix::t(SingleCellExperiment::logcounts(sds))) %>%
    dynwrap::add_trajectory(milestone_network = cluster_network,
                                progressions = progressions) %>% 
    dynwrap::add_dimred(dimred = SingleCellExperiment::reducedDims(sds)[[dimred_name]]) %>%
    dynwrap::add_grouping(grouping = grouping) %>%
    dynwrap::add_root(root_milestone_id = lineage_ctrl$start.clus) %>%
    dynwrap::add_pseudotime()
  output
}

data_prep_slingshot <- function(sobj, filename, reductions, path, cell_identities, start_cell){
  fs::path(path,"slingshot") %>% fs::dir_create()
  sce <- Seurat::as.SingleCellExperiment(sobj)
  start_cluster <- SingleCellExperiment::colData(sce)[start_cell,cell_identities] %>% as.character()
  available.reductions <- purrr::map_chr(reductions,~stringr::str_extract(names(SingleCellExperiment::reducedDims(sce)), .x) %>% purrr::discard(is.na))
  purrr::map(available.reductions, ~ list("clusterLabels" = cell_identities,
                                          "reducedDim" = .,
                                          "start.clus" = start_cluster,
                                          "end.clus" = NULL,
                                          "shrink" =  TRUE,
                                          "reweight" = TRUE,
                                          "reassign" = TRUE,
                                          "thresh" = 0.001,
                                          "maxit" = 10L,
                                          "stretch" = 2L,
                                          "approx_points" = 100L,
                                          "smoother" = "smooth.spline",
                                          "shrink.method" = "cosine"
                                          )) %>% 
    purrr::set_names(available.reductions) %>%
    yaml::write_yaml(fs::path(path,"slingshot","slingshot_params.yaml"))
  saveRDS(sce, file = fs::path(path,"slingshot",paste0(filename,'.rds'))) 
}

run_slingshot <- function(path, filename){
  path %<>% fs::path("slingshot")
  sce <- readRDS(fs::path(path,paste0(filename,'.rds')))
  slingshot_params <- yaml::read_yaml(fs::path(path,"slingshot_params.yaml")) %>%
    filter_computed_traj(path)
  
  create_trajectory_dir(slingshot_params,path)
  
  trajectory_sce <- purrr::map(slingshot_params, ~ list("traj" = slingshot::slingshot(data = sce, 
                                                                    clusterLabels = .[['clusterLabels']],
                                                                    reducedDim = .[['reducedDim']],
                                                                    start.clus = unlist(.[['start.clus']]), # to pass is.null()
                                                                    end.clus = unlist(.[['end.clus']]),
                                                                    shrink = .[['shrink']],
                                                                    reweight = .[['reweight']],
                                                                    reassign = .[['reassign']],
                                                                    thresh = .[['thresh']],
                                                                    maxit = .[['maxit']],
                                                                    stretch = .[['stretch']],
                                                                    approx_points = .[['approx_points']],
                                                                    smoother = .[['smoother']],
                                                                    shrink.method = .[['shrink.method']]
                                                                ),
                                                "redim" = .[['reducedDim']],
                                                "grouping" = .[['clusterLabels']]
                                                )
                       )  
  purrr::map(trajectory_sce,"traj") %>% purrr::iwalk( ~ saveRDS(slingshot::SlingshotDataSet(.x), file = fs::path(path,.y,"trajectory_result_slingshot_dataset.rds")))
  output <- purrr::map(trajectory_sce, ~ wrapper_sligshot(.[['traj']],stringr::str_to_upper(.[['redim']]), .[["grouping"]]))
  return(output)
}

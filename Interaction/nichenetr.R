run_nichnetr <- function(sobj, sender, receiver, condition_oi, condition_reference, network_name = "mouse", network_path, path, test = NULL, pval = 0.05, logFC = 0.58, expression_pct = 0.10, top_n_ligands = 20, top_n_targets = 200, cutoff_visualization = 0.33){
  require(ggplot2)
  assertthat::assert_that(is_seurat(sobj))
  assertthat::assert_that(is.character(sender), is.character(receiver))
  assertthat::assert_that(assertthat::is.dir(network_path), assertthat::is.readable(network_path))
  
  supported_network <- fs::dir_ls(network_path, type = "directory") %>% fs::path_file()
  assertthat::assert_that(is_in(network_name, supported_network))
  
  
  ## Network
  ligand_target_matrix <- fs::path(network_path, network_name, "ligand_target_matrix.rds") %>% readRDS()
  lr_network <- fs::path(network_path, network_name, "lr_network.rds") %>% readRDS()
  weighted_networks_lr <- fs::path(network_path, network_name, "weighted_networks_lr.rds") %>% readRDS()
  
  lr_network_strict <- lr_network %>% dplyr::filter(database != "ppi_prediction_go" & database != "ppi_prediction")
  
  ligands <- lr_network %>% dplyr::pull(from) %>% unique()
  receptors <- lr_network %>% dplyr::pull(to) %>% unique()
  ligands_bona_fide <- lr_network_strict %>% dplyr::pull(from) %>% unique()
  receptors_bona_fide <- lr_network_strict %>% dplyr::pull(to) %>% unique()
  
  expressed_genes_receiver <- receiver %>%
    unique() %>%
    purrr::set_names(.,.) %>%
    purrr::map(nichenetr::get_expressed_genes, seurat_obj = sobj, pct = 0.10) %>%
    unlist() %>%
    unique()
  background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  
  list_expressed_genes_sender <- sender %>%
    unique() %>%
    purrr::set_names(.,.) %>%
    purrr::map(nichenetr::get_expressed_genes, seurat_obj = sobj, pct = 0.10) 
    
  expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
  
  # DE table
  DE_type <- "conditions"
  
  available_test <- fs::path(path, "DIFFERENTIAL_EXPRESSION", DE_type) %>%
    fs::dir_ls() %>%
    fs::path_file()
  
  if (!is.null(test)) {
    assertthat::assert_that(is.character(test))
  } else {
    test <- available_test
  }
  
  if (any(!test %in% available_test)) {
    wrong_test <- test[!test %in% available_test]
    test <- setdiff(test, wrong_test)
    warning("The following test are not available, they will be ignored :\n", stringr::str_c(" - ", wrong_test, collapse = "\n"))
  }
  
  DE_table_receiver <- purrr::map_df(test, ~ fs::path(path, "DIFFERENTIAL_EXPRESSION", DE_type, .x, paste(.x, DE_type, "markers", sep = "_"), ext = "txt") %>%
                                       readr::read_tsv(col_types = test_cols_spec(.x)) %>%
                                       dplyr::select(gene, logFC, adj.P.Val, tested_cluster, control_cluster, cluster_conditions)) %>%
    dplyr::filter(tested_cluster == condition_oi, control_cluster == condition_reference, cluster_conditions %in% receiver, adj.P.Val < pval, abs(logFC) <= !!logFC) %>%
    dplyr::distinct(gene, .keep_all = TRUE)
  
  geneset_oi <- dplyr::pull(DE_table_receiver, gene)
  
  geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
  if (length(geneset_oi) == 0){
    stop("No genes were differentially expressed")
  }
  
  background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  
  # step3 nichenet analysis: define potential ligands
  expressed_ligands <- intersect(ligands,expressed_genes_sender)
  expressed_receptors <- intersect(receptors,expressed_genes_receiver)
  if (length(expressed_ligands) == 0){
    stop("No ligands expressed in sender cell")
  }
  if (length(expressed_receptors) == 0){
    stop("No receptors expressed in receiver cell")
  }
  potential_ligands <- lr_network %>% dplyr::filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
  if (length(potential_ligands) == 0){
    stop("No potentially active ligands")
  }
  
  ligand_activities <- nichenetr::predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
  ligand_activities %<>%
    dplyr::arrange(-pearson) %>%
    dplyr::mutate(rank = rank(desc(pearson)),
           bona_fide_ligand = test_ligand %in% ligands_bona_fide)
  
  if (top_n_ligands > 0){
    best_upstream_ligands <- ligand_activities %>% dplyr::top_n(top_n_ligands, pearson) %>% dplyr::arrange(-pearson) %>% dplyr::pull(test_ligand) %>% unique()
  }
  
  ligand_exp_plot <- Seurat::DotPlot(sobj, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + Seurat::RotatedAxis()
  
  active_ligand_target_links_df <- best_upstream_ligands %>% 
    purrr::map_dfr(nichenetr::get_weighted_ligand_target_links, geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = top_n_targets) %>%
    tidyr::drop_na()
  
  if(nrow(active_ligand_target_links_df) > 0){
    active_ligand_target_links <- nichenetr::prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = cutoff_visualization)
    order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
    order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
    rownames(active_ligand_target_links) <- rownames(active_ligand_target_links) %>% make.names()
    colnames(active_ligand_target_links) <- colnames(active_ligand_target_links) %>% make.names()
    
    order_targets %<>% intersect(rownames(active_ligand_target_links))
    order_ligands %<>% intersect(colnames(active_ligand_target_links))
    
    vis_ligand_target <- active_ligand_target_links[order_targets,order_ligands] %>% t()
    p_ligand_target_network <- vis_ligand_target %>% 
      nichenetr::make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) #+ scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.006,0.012))
  } else {
    vis_ligand_target = NULL
    p_ligand_target_network = NULL
    message("no highly likely active targets found for top ligands")
  }
  
  ligand_pearson_matrix <- ligand_activities %>% dplyr::select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
  
  rownames(ligand_pearson_matrix) <- rownames(ligand_pearson_matrix) %>% make.names()
  colnames(ligand_pearson_matrix) <- colnames(ligand_pearson_matrix) %>% make.names()
  
  vis_ligand_pearson <- ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
  p_ligand_pearson <- vis_ligand_pearson %>% nichenetr::make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9))
  
  combined_plot <- patchwork::wrap_plots(p_ligand_pearson, p_ligand_target_network + ylab(NULL), ncol = 2, widths = c(ncol(vis_ligand_pearson)+10, ncol(vis_ligand_target)), guides = "collect" ) &
    theme(legend.direction = "vertical")
  
  lr_network_top <- lr_network %>% dplyr::filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% dplyr::distinct(from,to)
  best_upstream_receptors <- lr_network_top %>% dplyr::pull(to) %>% unique()
  
  lr_network_top_df_large <- weighted_networks_lr %>% dplyr::filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
  
  lr_network_top_df <- lr_network_top_df_large %>% tidyr::pivot_wider(names_from = "from", values_from = "weight", values_fill = 0)
  lr_network_top_matrix <- lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
  
  if (nrow(lr_network_top_matrix) > 1){
    dist_receptors <- dist(lr_network_top_matrix, method = "binary")
    hclust_receptors <- hclust(dist_receptors, method = "ward.D2")
    order_receptors <- hclust_receptors$labels[hclust_receptors$order]
  } else {
    order_receptors <- rownames(lr_network_top_matrix)
  }
  if (ncol(lr_network_top_matrix) > 1) {
    dist_ligands <- dist(lr_network_top_matrix %>% t(), method = "binary")
    hclust_ligands <- hclust(dist_ligands, method = "ward.D2")
    order_ligands_receptor <- hclust_ligands$labels[hclust_ligands$order]
  } else {
    order_ligands_receptor <- colnames(lr_network_top_matrix)
  }
  
  order_receptors <- order_receptors %>% intersect(rownames(lr_network_top_matrix))
  order_ligands_receptor <- order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))
  
  vis_ligand_receptor_network <- lr_network_top_matrix[order_receptors, order_ligands_receptor]
  dim(vis_ligand_receptor_network) <- c(length(order_receptors), length(order_ligands_receptor))
  rownames(vis_ligand_receptor_network) <- order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network) <- order_ligands_receptor %>% make.names()
  
  p_ligand_receptor_network <- vis_ligand_receptor_network %>% t() %>% nichenetr::make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential") + theme(legend.position = "right")
  
  # bona fide ligand-receptor
  lr_network_top_df_large_strict <- lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
  lr_network_top_df_large_strict <- lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))
  
  lr_network_top_df_strict <- lr_network_top_df_large_strict %>% tidyr::pivot_wider(names_from = "from", values_from = "weight", values_fill = 0)
  lr_network_top_matrix_strict <- lr_network_top_df_strict %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)
  
  if (nrow(lr_network_top_df_large_strict) == 0){
    message("Remark: no bona fide receptors of top ligands")
    vis_ligand_receptor_network_strict <- NULL
    p_ligand_receptor_network_strict <- NULL
    lr_network_top_df_large_strict <-  NULL
    
  } else {
    if (nrow(lr_network_top_matrix_strict) > 1){
      dist_receptors <- dist(lr_network_top_matrix_strict, method = "binary")
      hclust_receptors <- hclust(dist_receptors, method = "ward.D2")
      order_receptors <- hclust_receptors$labels[hclust_receptors$order]
    } else {
      order_receptors <- rownames(lr_network_top_matrix)
    }
    if (ncol(lr_network_top_matrix_strict) > 1) {
      dist_ligands <- dist(lr_network_top_matrix_strict %>% t(), method = "binary")
      hclust_ligands <- hclust(dist_ligands, method = "ward.D2")
      order_ligands_receptor <- hclust_ligands$labels[hclust_ligands$order]
    } else {
      order_ligands_receptor <- colnames(lr_network_top_matrix_strict)
    }
    order_receptors <- order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
    order_ligands_receptor <- order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))
    
    vis_ligand_receptor_network_strict <- lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
    dim(vis_ligand_receptor_network_strict) <- c(length(order_receptors), length(order_ligands_receptor))
    
    rownames(vis_ligand_receptor_network_strict) <- order_receptors %>% make.names()
    colnames(vis_ligand_receptor_network_strict) <- order_ligands_receptor %>% make.names()
    
    p_ligand_receptor_network_strict <- vis_ligand_receptor_network_strict %>% t() %>% nichenetr::make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)") + theme(legend.position = "right")
    
    lr_network_top_df_large_strict <- lr_network_top_df_large_strict %>% dplyr::rename(ligand = from, receptor = to)
  }
  
  # list(p_ligand_receptor_network_strict,p_ligand_receptor_network, combined_plot)
  path %<>% fs::path("INTERACTION", "nichenetr", paste0("sender_", paste0(sender, collapse = "-"), "_receiver_",  paste0(receiver, collapse = "-"))) %>%
    fs::dir_create()
  # saving 
  ggsave(plot = combined_plot, filename = "regulatory_potential.png" , path = path, width = 20, height = 8)
  ggsave(plot = p_ligand_receptor_network, filename = "interaction_potential.png", path = path, width = 12, height = 9)
  ggsave(plot = p_ligand_receptor_network_strict, filename = "interaction_potential_strict.png", path = path , width = 8, height = 6)
  
  save(ligand_activities,best_upstream_ligands, active_ligand_target_links_df, lr_network_top_df_large,vis_ligand_target, active_ligand_target_links_df,
       vis_ligand_receptor_network, lr_network_top_df_large, vis_ligand_receptor_network_strict, lr_network_top_df_large_strict, geneset_oi, background_expressed_genes,
       file = "data_res_nichenetr.rda")
}

ligand_target_signaling_path <- function(ligands, targets, network_path, network_name){
  
  ligand_tf_matrix <- fs::path(network_path, network_name, "ligand_tf_matrix.rds") %>% readRDS()
  weighted_networks <- fs::path(network_path, network_name, "weighted_networks.rds") %>% readRDS()
  lr_network <- fs::path(network_path, network_name, "lr_network.rds") %>% readRDS()
  sig_network <- fs::path(network_path, network_name, "sig_network.rds") %>% readRDS()
  gr_network <- fs::path(network_path, network_name, "gr_network.rds") %>% readRDS()
  
  active_signaling_network <- nichenetr::get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, 
                                                                   ligands_all = ligands, 
                                                                   targets_all = targets, 
                                                                   weighted_networks = weighted_networks)
  
  active_signaling_network_scaled <- active_signaling_network %>% 
    purrr::map(~dplyr::mutate(.x, weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75))
  
  
  nichenetr::diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network_scaled, ligands_all = ligands, targets_all = targets, sig_color = "indianred", gr_color = "steelblue")
  
  data_source_network <- nichenetr::infer_supporting_datasources(signaling_graph_list = active_signaling_network,lr_network = lr_network, sig_network = sig_network, gr_network = gr_network)
  head(data_source_network) 
}

Seurat_test <- function(DE_type,output_path,test.use, sobj, name_first_group, name_second_group, min.pct, batch_list, verbose, cluster_conditions){
  #tests DE
  markers <- purrr::quietly(Seurat::FindMarkers)(sobj, ident.1 = name_first_group, ident.2 = name_second_group, min.pct = min.pct, logfc.threshold = 0, test.use = test.use, latent.vars = as.vector(batch_list), verbose = verbose)$result # MAST produce unwanted output
  colnames(markers) <- c("p_val","logFC","pct.1","pct.2","adj.P.Val")
  
  top20genes <- markers  %>%
    tibble::rownames_to_column(var = "gene") %>%
    dplyr::filter(adj.P.Val < 0.05) %>% 
    dplyr::top_n(20, abs(logFC)) %>% 
    dplyr::pull(gene)
  
  clusters <- Seurat::Idents(sobj) %>% tibble::enframe(name = "cells", value = "cluster")
  cells <- clusters %>% dplyr::filter(cluster == name_first_group | cluster == name_second_group) %>% dplyr::pull(cells)
  raw_counts <- Seurat::GetAssayData(object = sobj, slot = 'counts', assay = 'RNA')[top20genes,cells]#récup des expressions pour les top20genes
  log_count <- log(raw_counts+1)
  
  common_to_all_test(DE_type = DE_type, test = test.use, markers = markers, log_count = log_count, clusters = clusters, top20genes = top20genes, output_path = output_path, name_first_group = name_first_group, name_second_group = name_second_group, min.pct = min.pct, cluster_conditions = cluster_conditions)
}

Limma_voom_test <- function(DE_type,output_path, counts, clusters, design, contr.matrix, name_first_group, name_second_group, min.pct, cluster_conditions){
  #voom
  png(filename = fs::path(output_path, "Limma_Voom","Voom_plot",paste0("Limmavoom-",name_first_group,"-vs-",name_second_group,"-for",min.pct,"-plotVoom.png")))
  v <- limma::voom(counts, design = design, plot = T)
  dev.off()
  #exemple d'un plot valide et d'un plot non-valide là: https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
  #Limma Voom
  fitLV <-limma::lmFit(v, design)
  fitLimmaVoom <- limma::contrasts.fit(fitLV, contrasts=contr.matrix)
  fitLimmaVoom <-limma::eBayes(fitLimmaVoom, robust=FALSE)
  resLimmaVoom <-limma::topTable(fitLimmaVoom, number = nrow(counts),sort.by="logFC", adjust.method="BH")

  top20genes <- resLimmaVoom  %>%
    tibble::rownames_to_column(var = "gene") %>%
    dplyr::filter(adj.P.Val < 0.05) %>% 
    dplyr::top_n(20, abs(logFC)) %>% 
    dplyr::pull(gene)
  
  log_count <- log(counts$counts[as.character(top20genes),]+1)
  common_to_all_test(DE_type = DE_type, test = "Limma_Voom", markers = resLimmaVoom, log_count = log_count, clusters = clusters, top20genes = top20genes, output_path = output_path, name_first_group = name_first_group, name_second_group = name_second_group, min.pct = min.pct,cluster_conditions = cluster_conditions)
  
}

Limma_trend_test <- function(DE_type,output_path, counts, clusters, design, contr.matrix, name_first_group, name_second_group, min.pct, cluster_conditions){
  #Limma Trend
  logCPM <- edgeR::cpm(counts, log=TRUE, prior.count=0.5)# prior.count= ce qu'il faut rajouter au counts pour ne pas avoir log(0).
  fitLT <- limma::lmFit(logCPM, design=design)
  fitLimmaTrend <- limma::contrasts.fit(fitLT, contrasts=contr.matrix)
  fitLimmaTrend <- limma::eBayes(fitLimmaTrend, trend=TRUE, robust=FALSE)
  png(filename = fs::path(output_path, "Limma_Trend","Trend_plot",paste0("Limmatrend-",name_first_group,"-vs-",name_second_group,"-for",min.pct,"-plotTrend.png")))
  limma::plotSA(fitLimmaTrend, xlab = "Average log-expression", ylab = "sqrt(sigma)", zero.weights = FALSE)
  dev.off()
  resLimmaTrend <-limma::topTable(fitLimmaTrend, number = nrow(counts),sort.by="logFC", adjust.method="BH")
  
  top20genes <- resLimmaTrend  %>%
    tibble::rownames_to_column(var = "gene") %>%
    dplyr::filter(adj.P.Val < 0.05) %>% 
    dplyr::top_n(20, abs(logFC)) %>% 
    dplyr::pull(gene)
  
  log_count <- log(counts$counts[as.character(top20genes),]+1)
  common_to_all_test(DE_type = DE_type, test = "Limma_Trend", markers = resLimmaTrend, log_count = log_count, clusters = clusters, top20genes = top20genes, output_path = output_path, name_first_group = name_first_group, name_second_group = name_second_group, min.pct = min.pct, cluster_conditions = cluster_conditions)
}

# EdgeR_test <- function(DE_type, edge_type, output_path, counts, clusters, design, contr.matrix, name_first_group, name_second_group, min.pct){
#   output_path_plot <- switch(edge_type, 
#                              EdgeR_LRT = fs::path(output_path,"EdgeR_LRT"),
#                              EdgeR_QL = fs::path(output_path, "EdgeR_QL"))
#   counts <- switch(edge_type,
#                    EdgeR_LRT = edgeR::estimateGLMRobustDisp(counts, design), #long 40 min à 24Go
#                    EdgeR_QL = edgeR::estimateDisp(counts, design, robust=TRUE))
#   
#   png(filename = fs::path(output_path_plot, "BCV_plot", paste0(edge_type,"-",name_first_group,"-vs-",name_second_group,"-for",min.pct,"-plotBCV.png")))
#   edgeR::plotBCV(counts)
#   dev.off()
#   
#   if(edge_type == "EdgeR_LRT"){
#     resEdgeR <- edgeR::glmFit(counts, design, robust=TRUE) %>%
#       edgeR::glmLRT(contrast=contr.matrix) %>%
#       edgeR::topTags( adjust.method="BH", n=Inf, sort.by = 'none')
#   } else if (edge_type == "EdgeR_QL"){
#     fitEdgeRQL <- edgeR::glmQLFit(counts, design, robust=TRUE)
#     png(filename = fs::path(output_path_plot, "QLDisp_plot", paste0(edge_type,"-",name_first_group,"-vs-",name_second_group,"-for",min.pct,"-plotQLDisp.png")))
#     edgeR::plotQLDisp(fitEdgeRQL)
#     dev.off()
#     qlEdgeRQL<-edgeR::glmQLFTest(fitEdgeRQL, contrast=contr.matrix)
#     resEdgeR <- edgeR::topTags(qlEdgeRQL, adjust.method="BH", n=Inf, sort.by = 'none')
#   }
#   
#   top20genes <- resEdgeR  %>%
#     tibble::rownames_to_column(var = "gene") %>%
#     dplyr::filter(adj.P.Val < 0.05) %>% 
#     dplyr::top_n(20, abs(logFC)) %>% 
#     dplyr::pull(gene)
#   
#   
#   colnames(resEdgeR$table)[4]=c("adj.P.Val")
#   log_count <- log(counts$counts[as.character(top20genes),]+1)
# 
#   common_to_all_test(DE_type = DE_type, test = edge_type, markers = resEdgeR, log_count = log_count, clusters = clusters, top20genes = top20genes, output_path = output_path, name_first_group = name_first_group, name_second_group = name_second_group, min.pct = min.pct)
# }

EdgeR_LRT_test <- function(DE_type,output_path, counts, clusters, design, contr.matrix, name_first_group, name_second_group, min.pct, cluster_conditions){
  #edgeRLRT
  counts <- edgeR::estimateGLMRobustDisp(counts, design) #long 40 min à 24Go
  png(filename = fs::path(output_path,"EdgeR_LRT", "BCV_plot", paste0("EdgeR_LRT-",name_first_group,"-vs-",name_second_group,"-for",min.pct,"-plotBCV.png")))
  plotBCV(counts)
  dev.off()
  fitEdgeRLRT <- edgeR::glmFit(counts, design, robust=TRUE)
  lrtEdgeRLRT<-edgeR::glmLRT(fitEdgeRLRT, contrast=contr.matrix)
  resEdgeRLRT <- edgeR::topTags(lrtEdgeRLRT, adjust.method="BH", n=Inf, sort.by = 'none')
  
  colnames(resEdgeRLRT$table)[4]=c("adj.P.Val")
  
  top20genes <- resEdgeRLRT$table  %>%
    tibble::rownames_to_column(var = "gene") %>%
    dplyr::filter(adj.P.Val < 0.05) %>% 
    dplyr::top_n(20, abs(logFC)) %>% 
    dplyr::pull(gene)
  
  log_count <- log(counts$counts[as.character(top20genes),]+1)
  
  common_to_all_test(DE_type = DE_type, test = "EdgeR_LRT", markers = resEdgeRLRT$table, log_count = log_count, clusters = clusters, top20genes = top20genes, output_path = output_path, name_first_group = name_first_group, name_second_group = name_second_group, min.pct = min.pct, cluster_conditions = cluster_conditions)
}

EdgeR_QL_test <- function(DE_type,output_path, counts, clusters, design, contr.matrix, name_first_group, name_second_group, min.pct, cluster_conditions){
  #edgeRQL
  counts <- edgeR::estimateDisp(counts, design, robust=TRUE)
  png(filename = fs::path(output_path,"EdgeR_QL", "BCV_plot", paste0("EdgeR_QL-",name_first_group,"-vs-",name_second_group,"-for",min.pct,"-plotBCV.png"))) 
  edgeR::plotBCV(counts)
  dev.off()
  fitEdgeRQL <- edgeR::glmQLFit(counts, design, robust=TRUE)
  png(filename = fs::path(output_path, "EdgeR_QL", "QLDisp_plot", paste0("EdgeR_QL-",name_first_group,"-vs-",name_second_group,"-for",min.pct,"-plotQLDisp.png")))
  edgeR::plotQLDisp(fitEdgeRQL)
  dev.off()
  qlEdgeRQL<-edgeR::glmQLFTest(fitEdgeRQL, contrast=contr.matrix)
  resEdgeRQL <- edgeR::topTags(qlEdgeRQL, adjust.method="BH", n=Inf, sort.by = 'none')
  
  colnames(resEdgeRQL$table)[4] <- c("adj.P.Val")
  
  top20genes <- resEdgeRQL$table  %>%
    tibble::rownames_to_column(var = "gene") %>%
    dplyr::filter(adj.P.Val < 0.05) %>% 
    dplyr::top_n(20, abs(logFC)) %>% 
    dplyr::pull(gene)
  
  
  
  log_count <- log(counts$counts[as.character(top20genes),]+1)
  
  common_to_all_test(DE_type = DE_type, test = "EdgeR_QL", markers = resEdgeRQL$table, log_count = log_count, clusters = clusters, top20genes = top20genes, output_path = output_path, name_first_group = name_first_group, name_second_group = name_second_group, min.pct = min.pct, cluster_conditions = cluster_conditions)
}

common_to_all_test <- function(DE_type, test, markers, log_count, clusters, top20genes, output_path, name_first_group, name_second_group, min.pct, cluster_conditions){
  if(DE_type == "1vsAll"){name_second_group2  <- "All"} else {name_second_group2  <- name_second_group}
  # annot <- data.frame(sample=clusters[colnames(log_count),]) 
  # rownames(annot) <- colnames(log_count)
  
  # Display a warning or comment on the plot and check top markers
  # any(is.infinite(markers[['logFC']])) 
  filename <- format_filename(min.pct, name_first_group, name_second_group2, cluster_conditions)
  output_path %<>% fs::path(test)
  markers %<>% tibble::rownames_to_column(var = "gene") %>% dplyr::filter(is.finite(logFC)) %>%
    tibble::add_column(tested_cluster = name_first_group, control_cluster = name_second_group, min.pct = min.pct, cluster_conditions = stringr::str_c(cluster_conditions, collapse = "-"))
  write_DE_test(df = markers, filename = file.path(output_path,paste0(test,"_",DE_type,"_markers.txt")))
  
  
  plot_pval <- ggplot2::ggplot(markers, ggplot2::aes(x=adj.P.Val)) + ggplot2::geom_histogram(bins = 30) + ggplot2::labs(title = stringr::str_replace_all(test,"_"," "), x="Adjusted P-values")
  suppressMessages(ggplot2::ggsave(filename = paste0(filename,"-plotPVal.png"), plot = plot_pval, path = fs::path(output_path, "pvalues_plot")))
  
  if (!rlang::is_empty(top20genes)) {
    #heatmap
    plot_heatmap <- heatmap(log_count = log_count, clusters = clusters, DE_type = DE_type, row_tree = FALSE, annotation = TRUE, label = FALSE, name_first_group = name_first_group, name_second_group = name_second_group2) + patchwork::plot_annotation(title = stringr::str_replace_all(test,"_"," "), subtitle = paste0(name_first_group," vs ",name_second_group2," for min.pct : ",min.pct*100," %"))
    suppressMessages(ggplot2::ggsave(filename = paste0(filename,"-plotHeatmap.png"), plot = plot_heatmap, path = fs::path(output_path,"heatmap_plot")))
  }
 
  #volcanoplot
  logFC_amp <- markers[["logFC"]] %>%
    summary() %>%
    .[c("Min.", "Max.")] %>%
    diff()
  
  if (logFC_amp > 10) {
    label_text <- ggplot2::geom_text(mapping = ggplot2::aes(x = Inf, y = -log10(0.05), label = "P-value cut off 0.05", vjust = 1.5, hjust = 1))
  } else {
    label_text <- list(
      ggplot2::geom_text(mapping = ggplot2::aes(x = 0.58, y = Inf, label = "0.58", vjust = -0.5)),
      ggplot2::geom_text(mapping = ggplot2::aes(x = -0.58, y = Inf, label = "- 0.58", vjust = -0.5)),
      ggplot2::geom_text(mapping = ggplot2::aes(x = Inf, y = -log10(0.05), label = "P-value cut off 0.05", vjust = 1.5, hjust = 1))
    )
  }
  
  plot_volcano <- EnhancedVolcano::EnhancedVolcano(markers, lab = markers$gene, x = "logFC", y = "adj.P.Val", axisLabSize = 14, title = stringr::str_replace_all(test, "_", " "), subtitle = paste0(name_first_group, " vs ", name_second_group2, " for min.pct : ", min.pct * 100, " %"), ylab = bquote(~ -Log[10] ~ (adjusted ~ P - values)), legendPosition = "bottom", selectLab = top20genes, pCutoff = 0.05, FCcutoff = 0.58, gridlines.major = FALSE, gridlines.minor = FALSE, boxedLabels = TRUE, drawConnectors = TRUE, ylim = c(-50, 400)) + 
    label_text + 
    ggplot2::coord_cartesian(clip = "off")
  suppressMessages(ggplot2::ggsave(paste0(filename,"-plotVolcano.png"), plot = plot_volcano, path = fs::path(output_path,"volcano_plot")))
  
}

calc_features_pct <- function(list_cell, min.pct) {
  thresh.min <- 0 #si on veut au moins 1 comptage par cellule, on met "thresh.min <- 0"
  data.alpha=NULL
  
  for(i in 1:length(list_cell)){
    #CALCULS DES POURCENTAGES EXPRESSION POUR CHQ GENE DE CHQ GRP
    pct <- round(Matrix::rowSums(x = (list_cell[[i]] > thresh.min))/length(x = colnames(list_cell[[i]])), digits = 3) #attention la précison est celle implémentée dans Seurat, digit=3 correspond à 0,000 (donc 0,0%)
    if (is.null(data.alpha)){
      data.alpha=pct
    }else{
      data.alpha=cbind(data.alpha,pct)#tableau avec les pourcentages
    }
  }
  colnames(x = data.alpha) <- c("pct.1", "pct.2")
  #CALCULER LES GENES QUI PASSENT LE THRESHOLD
  alpha.min <- apply(X = data.alpha, MARGIN = 1, FUN = max)#tableau avec le pourcentage le plus élevé entre les groupes comparés
  names(x = alpha.min) <- rownames(x = data.alpha)
  #features   <- names(x = which(x = alpha.min > min.pct))
  features=data.alpha[which(x = alpha.min > min.pct),]
  features
}

#attribut de nouveau nom de cluster à un objet seurat en fonction des groupes de cluster à comparer (nécessaire pour les tests Seurat puisqu'ils n'utilisent pas de contraste)
calc_new.cluster <- function(sobj, first_group, second_group, list_element){
  new.cluster=NULL
  for (element in list_element){
    if(element %in% first_group){
      if (is.null(new.cluster)){
        new.cluster=str_c(first_group,collapse = "_") #name_first_group
      }else{
        new.cluster=c(new.cluster,str_c(first_group,collapse = "_"))
      }
    }else if (element %in% second_group){
      if (is.null(new.cluster)){
        new.cluster=str_c(second_group,collapse = "_") #name_second_group
      }else{
        new.cluster=c(new.cluster,str_c(second_group,collapse = "_"))
      }
    }else{ #on ne change pas le nom du cluster
      if (is.null(new.cluster)){
        new.cluster=element
      }else{
        new.cluster=c(new.cluster,element)
      }
    }
  }
  names(new.cluster) <- levels(sobj)
  sobj <- RenameIdents(sobj, new.cluster)
  sobj
}

#Calcul les gènes restant après le min.pct dans les 2 groupes à comparer
# calc_features_pct <- function(counts, cells.1, cells.2, min.pct)  {
#   #CALCULS DES POURCENTAGES EXPRESSION POUR CHQ GENE DE CHQ GRP
#   thresh.min <- 0 #si on veut au moins 1 comptage par cellule, on met "thresh.min <- 0"
#   pct.1 <- round(Matrix::rowSums(x = (cells.1 > thresh.min))/length(x = colnames(cells.1)), digits = 3) #attention la précison est celle implémentée dans Seurat, digit=3 correspond à 0,000 (donc 0,0%)
#   pct.2 <- round(Matrix::rowSums(x = (cells.2 > thresh.min))/length(x = colnames(cells.2)), digits = 3)
#   #CALCULER CE QUI PASSE LE THRESHOLD
#   data.alpha <- cbind(pct.1, pct.2) #tableau avec les pourcentages
#   colnames(x = data.alpha) <- c("pct.1", "pct.2")
#   alpha.min <- apply(X = data.alpha, MARGIN = 1, FUN = max)#tableau avec le pourcentage le plus élevé entre les deux conditions
#   names(x = alpha.min) <- rownames(x = data.alpha)
#   features   <- names(x = which(x = alpha.min > min.pct))
#   features
# }

#attribut de nouveau nom de cluster à un objet seurat en fonction des groupes de cluster à comparer (nécessaire pour les tests Seurat puisqu'ils n'utilisent pas de contraste)
# calc_new.cluster <- function(sobj, first_group, second_group){
#   new.cluster <- NULL
#   clusters_list <- sort(unique(as.data.frame(Seurat::Idents(sobj))[,1]))
#   if (length(first_group) != 1 | length(second_group) != 1){
#     for (cluster in clusters_list){
#       if(cluster %in% first_group){
#         if (is.null(new.cluster)){
#           new.cluster <- stringr::str_c(first_group,collapse = "_") #name_first_group
#         }else{
#           new.cluster <- c(new.cluster,stringr::str_c(first_group,collapse = "_"))
#         }
#       }else{
#         if (is.null(new.cluster)){
#           new.cluster <- stringr::str_c(second_group,collapse = "_") #name_second_group
#         }else{
#           new.cluster <- c(new.cluster,stringr::str_c(second_group,collapse = "_"))
#         }
#       }
#     }
#     names(new.cluster) <- levels(sobj)
#     sobj <- Seurat::RenameIdents(sobj, new.cluster)
#   }
#   sobj
# }

#fabrication de la matrice de contraste
calc_contr.matrix <- function(first_group, second_group, design)  {
  c <- NULL
  for(i in first_group){
    if (is.null(c)){
      c <- paste0("C_",i)
    }else{
      c <- paste0(c,"+C_",i)
    }
  }
  d <- NULL
  for (i in second_group){
    if (is.null(d)){
      d <- paste0("C_",i)
    }else{
      d <- paste0(d,"+C_",i)
    }
  }
  contrast <- paste0("((",c, ")/",length(first_group),")-((", d, ")/", length(second_group), ")")
  contr.matrix <- limma::makeContrasts(contrasts=contrast,levels=design)
  contr.matrix
}

create_folders_DE <- function(out.dir, test){
  plots <- c('pvalues_plot','heatmap_plot','volcano_plot','BCV_plot','QLDisp_plot','Trend_plot','Voom_plot')
  
  purrr::map(plots,~fs::path(out.dir,test, .)) %>% unlist() %>%
    stringr::str_subset("Trend/Trend|Voom/Voom|QL/BCV|QL/QL|LRT/BCV|pvalues|heatmap|volcano") %>% 
    fs::dir_exists() %>% 
    .[. == FALSE] %>% 
    names() %>%
    fs::dir_create(recurse = TRUE)
  # purrr::walk(paths,function(p) if (!file.exists(p)){ dir.create(p,showWarnings = FALSE, recursive = TRUE) })
}

write_DE_test <- function(df, filename){
  file_lock <- paste0(filename,".lock")
  lock <- filelock::lock(path = file_lock, exclusive = TRUE, timeout = Inf)
  if (!fs::file_exists(filename)){
    readr::write_tsv(df, path = filename, quote_escape = FALSE, na = "NA", col_names = TRUE, append = FALSE)
  }else{
    readr::write_tsv(df, path = filename, quote_escape = FALSE, na = "NA", col_names = FALSE, append = TRUE)
  }
  filelock::unlock(lock)
}

generate_text_label <- function(DE_type, name_first_group, name_second_group){
  is_num <- stringr::str_split(name_first_group,"_") %>% unlist() %>% stringr::str_detect("^[0-9]+$") %>% all()
  if (DE_type == "1vsAll"){
    if(is_num){
      name <- stringr::str_c("Cluster : ", name_first_group)
    } else {
      name <- name_first_group
    }
    name <- c(name, "All")
  } else if (DE_type == "1vs1"){
    if(is_num){
      name <- stringr::str_c("Cluster : ", c(name_first_group, name_second_group))
    } else {
      name <- c(name_first_group, name_second_group)
    }
  } else if (DE_type == "SvsS"){
    if(is_num){
      name <- stringr::str_c("Cluster ", stringr::str_replace_all(c(name_first_group, name_second_group), "_", " - "))
    } else {
      name <- stringr::str_replace_all(c(name_first_group, name_second_group), "_", " - ")
    }
  } else if (DE_type == "conditions"){
    name <- stringr::str_replace_all(c(name_first_group, name_second_group), "_", " - ")
  }
  purrr::set_names(name,c(name_first_group, name_second_group)) %>% 
    tibble::enframe(name = 'cluster', value = 'name')
}


heatmap <- function(log_count, clusters, DE_type, row_tree = FALSE, annotation = TRUE, label = FALSE, name_first_group, name_second_group, colors = c("blue","yellow","red")){
  n_cells <- ncol(log_count)
  clust_levels <- levels(clusters$cluster)
  current_clust <- c(name_first_group, name_second_group)
  n_cluster <- nlevels(clusters$cluster)
  if(DE_type == "1vsAll" & name_second_group == "All"){
    clusters <- dplyr::mutate(clusters, cluster = ifelse(cluster == name_first_group, name_first_group, "All") ) 
  } else if (DE_type == "SvsS"){
    group_1_split <- stringr::str_split(name_first_group, "_") %>% unlist()
    group_2_split <- stringr::str_split(name_second_group, "_") %>% unlist()
    clusters <- dplyr::mutate(clusters, cluster = dplyr::case_when( cluster %in% group_1_split ~ name_first_group,
                                                                    cluster %in% group_2_split ~ name_second_group,
                                                                    TRUE ~ as.character(cluster)) ) %>%
      dplyr::mutate(cluster = as.factor(cluster))
  }

  clust_tree <- log_count %>% dist() %>% hclust() %>% as.dendrogram()
  tree <- ggdendro::dendro_data(clust_tree, type = "rectangle")
  ordered_gene <- ggdendro::label(tree) %>% dplyr::pull(label) %>% as.character()
  
  data <- tibble::as_tibble(log_count, rownames = "genes") %>%
    tidyr::pivot_longer(cols = -genes, names_to = "cells") %>% 
    dplyr::left_join(clusters, by = 'cells') %>%
    dplyr::arrange(cluster) %>%
    dplyr::mutate(cells = factor(cells) %>% forcats::fct_inorder(), genes = factor(genes) %>% forcats::fct_relevel(ordered_gene))
  
  
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = cells, y = genes, fill = value)) + 
    ggplot2::geom_tile() + 
    ggplot2::scale_fill_gradientn(name = "Log Expression", colours = colors, na.value = "white") +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(), 
                   axis.line = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank()) + 
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::theme(plot.margin = ggplot2::margin())
  
  genes <- levels(data$genes)
  n_genes <- nlevels(data$genes)
  
  axis <- ggplot2::ggplot(
    data.frame(
      y = 1:n_genes,
      gene = genes
    ),
    ggplot2::aes(x = 0, y = y, label = gene)
  ) +
    ggrepel::geom_text_repel(
      min.segment.length = grid::unit(0, "pt"),
      color = "grey30", ## ggplot2 theme_grey() axis text
      size = 0.8 * 11 / ggplot2::.pt ## ggplot2 theme_grey() axis text
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 0.5), expand = c(0, 0),
      breaks = NULL, labels = NULL, name = NULL
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0.5, n_genes + .5), expand = c(0, 0),
      breaks = NULL, labels = NULL, name = NULL
    ) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      plot.margin = ggplot2::margin()
    )
  
    # Preserve colors from seurat when 1 vs 1
    colors_annot <- scales::hue_pal()(n_cluster) %>% 
      magrittr::set_names(clust_levels) %>% 
      tibble::enframe(name = 'cluster', value = 'cols') %>% 
      dplyr::filter(cluster %in% current_clust) %>%
      as.data.frame()
    
    if(DE_type == "1vsAll"){
      if(nrow(colors_annot) == 1){
        colors_annot <- tibble::add_row(colors_annot, cluster = "All", cols = "#696969")
      } else {
        colors_annot[[2,1]] <- "All"
        colors_annot[[2,2]] <- "#696969"
      }
    } else if (DE_type == "SvsS") {
      colors_annot <- tibble::add_row(colors_annot, cluster = c(name_first_group, name_second_group), cols = c("#e41a1c","#377eb8"))
    } else if (DE_type == "conditions") {
      colors_annot <- tibble::add_row(colors_annot, cluster = c(name_first_group, name_second_group), cols = c("#e41a1c","#377eb8"))
    }
    
    # Create a color bar to represent the cells cluster
   color_test <- dplyr::group_by(clusters, cluster) %>% 
      dplyr::count() %>% 
      dplyr::filter(cluster == name_first_group | cluster == name_second_group) %>% 
      dplyr::ungroup(cluster) %>% 
      # dplyr::mutate(n = n /n_cells) %>%
      dplyr::inner_join(colors_annot, by = 'cluster') %>% 
      {c(rep(.[[1,3]],.[[1,2]]), rep(.[[2,3]],.[[2,2]]))} %>% t()
   
   bar <- ggplot2::ggplot() +
     ggplot2::theme_void() +
     ggplot2::annotation_raster(color_test, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
     ggplot2::scale_x_continuous(expand = c(0, 0)) +
     ggplot2::theme(plot.margin = ggplot2::margin())

    name <- generate_text_label(DE_type, name_first_group, name_second_group)
    
    txt1 <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = name[[1,2]], fontface = "bold") +
      ggplot2::scale_x_continuous(
        limits = c(0, 1), expand = c(0, 0),
        breaks = NULL, labels = NULL, name = NULL
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, 1), expand = c(0, 0),
        breaks = NULL, labels = NULL, name = NULL
      ) +
      ggplot2::theme(
        panel.background = ggplot2::element_blank(),
        plot.margin = ggplot2::margin()
      )
    
    txt2 <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label =  name[[2,2]], fontface = "bold") +
      ggplot2::scale_x_continuous(
        limits = c(0, 1), expand = c(0, 0),
        breaks = NULL, labels = NULL, name = NULL
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, 1), expand = c(0, 0),
        breaks = NULL, labels = NULL, name = NULL
      ) +
      ggplot2::theme(
        panel.background = ggplot2::element_blank(),
        plot.margin = ggplot2::margin()
      )
  
    dendro_row <- cowplot::axis_canvas(plot, axis = "y", coord_flip = TRUE, ylim = pbuild$layout$panel_params[[1]]$y.range) + 
      ggplot2::geom_segment(data = ggdendro::segment(tree), ggplot2::aes(y = -y, x = x, xend = xend, yend = -yend)) + 
      ggplot2::coord_flip() 
    
    patchwork::wrap_plots(txt1,
                          bar,
                          txt2,
                          dendro_row, 
                          plot, 
                          axis,
                          nrow = 2,
                          ncol = 3,
                          guides = 'collect',
                          heights = c(0.5,5), 
                          widths = c(1,5))
}

format_filename <- function(min.pct, name_first_group, name_second_group, cluster_conditions){
  name <- paste0(name_first_group,"-vs-",name_second_group,"-for-",min.pct)
  if (!is.null(cluster_conditions)){
    name <- paste0("clusters_",paste0(cluster_conditions, collapse = "-"),"-comparison_",name)
  } 
  name
}
stat_DE <- function(path, pval = 0.05, thresholdLogFC = 0.58, thresholdpct = 5, species = "mouse", verbose = TRUE){
  assertthat::assert_that(assertthat::is.dir(path), assertthat::is.writeable(path))
  assertthat::assert_that(is.numeric(pval), pval >= 0, pval <=1)
  assertthat::assert_that(is.numeric(thresholdpct), thresholdpct >= 0, thresholdpct <=100)
  assertthat::assert_that(is.numeric(thresholdLogFC))
  supported_species <- get_db()
  assertthat::assert_that(is.character(species), is_in(species, supported_species))
  
  path <- fs::path_tidy(path)
  output_path <- fs::path(path,"DE_STAT", paste0("pval", pval, "_thresholdLogFC",thresholdLogFC, "_thresholdpct", thresholdpct))
  
  DE_type <- fs::dir_ls(fs::path(path, "DIFFERENTIAL_EXPRESSION"), type = "directory") %>% fs::path_file() #c("1vsAll","1vs1","SvsS")
  
  if(verbose) message("Running DE statistics : ...")
  purrr::walk(DE_type,stat_by_vs, output_path = output_path, path = path, species = species, verbose = verbose, pval = pval, thresholdLogFC = thresholdLogFC, thresholdpct = thresholdpct)
}

stat_by_vs <- function(DE_type,output_path, path, species, verbose, pval, thresholdLogFC, thresholdpct){
  if(verbose) message("\t",DE_type)
  # create output folders
  output_path_Comparisons <- fs::path(output_path, DE_type, "DE_Comparisons")
  
  fs::path(output_path_Comparisons,c("nb_DE_genes","nb_common_DE_genes")) %>% 
    fs::dir_exists() %>% 
    .[. == FALSE] %>% 
    names() %>%
    fs::dir_create(recurse = TRUE)
  
  # available test for this DE_type
  test.used <- fs::dir_ls(fs::path(path, "DIFFERENTIAL_EXPRESSION", DE_type), type = "directory", recurse = FALSE) %>% fs::path_file()
  # Run for each test
  total_res <- purrr::map(test.used,stat_by_test, path = path, DE_type = DE_type, species = species, verbose = verbose, output_path = output_path, pval = pval, thresholdLogFC = thresholdLogFC, thresholdpct = thresholdpct)
  
  # Bind result together
  total_recap <- purrr::map_df(total_res,'recap')
  recap_DE_genes <- purrr::map_df(total_res,'recap_DE_genes') %>%
    dplyr::rowwise() %>%
    dplyr::mutate(control_cluster = ifelse(DE_type == "1vsAll", "All", as.character(control_cluster))) %>%
    dplyr::ungroup()
  rm(total_res)
  
  total_recap %<>% 
    dplyr::rowwise() %>%
    dplyr::mutate(control_cluster = ifelse(DE_type == "1vsAll", "All", as.character(control_cluster))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(x_ticks = stringr::str_c(tested_cluster,'vs',control_cluster, sep = " "))
  #plots
  n_comp <- total_recap %>% dplyr::select(tested_cluster,control_cluster, min.pct) %>% dplyr::distinct() %>% nrow()
  n_comp <- ifelse(n_comp < 20, 20L, n_comp)
  p <-  ggplot2::ggplot(total_recap,  ggplot2::aes(x = x_ticks, y=nb_sig_genes, fill=test)) +
    ggplot2::geom_dotplot(binaxis='y', stackratio=0.5, dotsize=0.4) +
    ggplot2::labs(title="Comparison of number of significatives genes", x="", y="Nb DE genes") +
    ggplot2::scale_fill_discrete(name = "") + 
    ggplot2::facet_wrap(~min.pct) + 
    ggplot2::theme(axis.text.x =  ggplot2::element_text(angle=90, vjust = 0.5))
  suppressMessages(ggplot2::ggsave(filename = "Comparison_nb_significatives_genes.png", plot = p, path = fs::path(output_path_Comparisons,"nb_DE_genes"), width = 0.3*n_comp + 1.4))
  
  p <-  ggplot2::ggplot(total_recap,  ggplot2::aes(x=x_ticks, y=nb_DE_genes, fill=test)) +
    ggplot2::geom_dotplot(binaxis='y', stackratio=0.5, dotsize=0.4) +
    ggplot2::labs(title="Comparison of number of DE genes between all DE tests", x="", y="Nb DE genes")+
    ggplot2::scale_fill_discrete(name = "") + 
    ggplot2::facet_wrap(~min.pct) + 
    ggplot2::theme(axis.text.x =  ggplot2::element_text(angle=90, vjust = 0.5))
  suppressMessages(ggplot2::ggsave(filename = "Comparison_nb_DEgenes.png", plot = p, path = fs::path(output_path_Comparisons,"nb_DE_genes"), width = 0.3*n_comp + 1.4))
  
    ##Comparaison des noms de gènes DE:
  contrast <- total_recap %>% dplyr::select("tested_cluster","control_cluster", "min.pct") %>% dplyr::distinct()
  
  if(verbose) {
    message("\tUpSet plotting")
    pb_upset <- progress::progress_bar$new(format = "\t\tPlotting : |:bar| :percent (:elapsedfull)",total = nrow(contrast), complete = "█", incomplete = " ", current = " ", clear = FALSE)
  } else {
    pb_upset <- NULL
  }
  
  purrr::pwalk(contrast, ~ recap_DE_genes %>% dplyr::filter(tested_cluster == ..1 & control_cluster == ..2 & min.pct == ..3) %>% plot_upset(tested_cluster = ..1, control_cluster = ..2, min.pct = ..3,test = test.used, output_path_Comparisons = output_path_Comparisons, pb = pb_upset))
  
}

plot_upset <- function(recap_DE_genes,tested_cluster, control_cluster, min.pct, test,  output_path_Comparisons, pb = NULL){
  upset_plot <- purrr::map(test,~dplyr::filter(recap_DE_genes, test_used == .) %>% 
                             dplyr::pull(genes) ) %>% purrr::set_names(test) %>% purrr::compact() %>% UpSetR::fromList() %>%
    UpSetR::upset(nsets = length(.), order.by = 'freq')
  png(file = fs::path(output_path_Comparisons,"nb_common_DE_genes",paste0("Comparison_names_DEgenes_",tested_cluster,"vs",control_cluster,"for",min.pct,".png")), width = 7, height = 5, res = 300, units = 'in')
  print(upset_plot)
  dev.off()
  if(!is.null(pb)) pb$tick()
}

stat_by_test <- function(test,path,DE_type, species, verbose, output_path, pval, thresholdLogFC, thresholdpct){
  # create folders
  output_path_def <- fs::path(output_path, DE_type, test)
  if (! fs::dir_exists(output_path_def)) fs::dir_create(output_path_def) 
  
  if(verbose) message("\tTest : ", test)
  
  # READ DATA
  if(any(c("wilcox","bimod","t","poisson","negbinom","LR","MAST", "roc", "DEseq2") == test)){
    data <- readr::read_tsv(file = fs::path(path,'DIFFERENTIAL_EXPRESSION',DE_type,test,paste0(test,"_",DE_type,"_markers.txt")), na = c("", "NA", "NaN"), col_names = TRUE, col_types = readr::cols(
      p_val = readr::col_number(),
      logFC = readr::col_number(),
      pct.1 = readr::col_number(),
      pct.2 = readr::col_number(),
      adj.P.Val = readr::col_number(),
      tested_cluster = readr::col_factor(),
      control_cluster = readr::col_factor(),
      gene = readr::col_character(),
      min.pct = readr::col_number(), 
      cluster_conditions = readr::col_character()
    ))
  } else if (stringr::str_detect(test,"EdgeR")){
    data <- readr::read_tsv(file = fs::path(path,'DIFFERENTIAL_EXPRESSION',DE_type,test,paste0(test,"_",DE_type,"_markers.txt")), na = c("", "NA", "NaN"), col_names = TRUE, col_types = readr::cols(
      logFC = readr::col_number(),
      logCPM = readr::col_number(),
      F = readr::col_number(),
      adj.P.Val = readr::col_number(),
      FDR = readr::col_number(),
      tested_cluster = readr::col_factor(),
      control_cluster = readr::col_factor(),
      gene = readr::col_character(),
      min.pct = readr::col_number(), 
      cluster_conditions = readr::col_character()
    ))
  } else if (stringr::str_detect(test,"Limma")){
    data <- readr::read_tsv(file = fs::path(path,'DIFFERENTIAL_EXPRESSION',DE_type,test,paste0(test,"_",DE_type,"_markers.txt")), na = c("", "NA", "NaN"), col_names = TRUE, col_types = readr::cols(
      logFC = readr::col_number(),
      AveExpr = readr::col_number(),
      t = readr::col_number(),
      P.Value = readr::col_number(),
      adj.P.Val = readr::col_number(),
      B = readr::col_number(),
      tested_cluster = readr::col_factor(),
      control_cluster = readr::col_factor(),
      gene = readr::col_character(),
      min.pct = readr::col_number(), 
      cluster_conditions = readr::col_character()
    ))
  }
  
  
  contrast <- data %>% dplyr::select("tested_cluster","control_cluster", "min.pct") %>% dplyr::distinct()
  
  if(verbose){
    pb_contrast <- progress::progress_bar$new(format = "\t\tComputing : |:bar| :percent (:elapsedfull)",total = nrow(contrast), complete = "█", incomplete = " ", current = " ", clear = FALSE)
  } else {
    pb_contrast <- NULL
  }
  
  res_by_contrast <- purrr::pmap(contrast,stat_by_contrast, data = data, species = species, pval = pval, thresholdLogFC = thresholdLogFC, thresholdpct = thresholdpct, pb = pb_contrast)
  
  recap <- purrr::map_df(res_by_contrast,"recap") %>% dplyr::mutate(pct_sig_not_found = (nb_sig_genes_not_found / nb_sig_genes) * 100,
                                                                    pct_DE_not_found = (nb_DE_genes_not_found / nb_DE_genes) * 100) 
  
  recap_DE_genes <- purrr::map_df(res_by_contrast,"recap_DE_genes") %>% tibble::add_column(test_used = test)
  recap <- recap %>% tibble::add_column(test = test)
  All_geneList_sig <- purrr::map(res_by_contrast,"All_geneList_sig")
  All_geneDE <- purrr::map(res_by_contrast,"All_geneDE")
  All_geneList <- purrr::map(res_by_contrast,"All_geneList")
  params <- purrr::map(res_by_contrast,"params")
  
  # Saving results
  readr::write_tsv(recap, path = fs::path(output_path_def, paste0(test,"_recap.txt")))
  readr::write_tsv(recap_DE_genes, path = fs::path(output_path_def, paste0(test,"_recap_DE_genes.txt")))
  saveRDS(All_geneDE, file = fs::path(output_path_def, paste0(test,"_all_geneDE.rds")))
  saveRDS(All_geneList, file = fs::path(output_path_def, paste0(test,"_All_geneList.rds")))
  
  # ENsure that fators are sorted in a numeric order
  recap %<>% dplyr::mutate(control_cluster = forcats::fct_relevel(control_cluster,stringr::str_sort(levels(control_cluster), numeric = TRUE)),
                           tested_cluster = forcats::fct_relevel(tested_cluster,stringr::str_sort(levels(tested_cluster), numeric = TRUE)),
                           color = factor(dplyr::if_else(pct_sig_not_found>thresholdpct | pct_DE_not_found>thresholdpct, "red", "black"),levels = c("black", "red")))
  n_more_percent <- recap %>% dplyr::group_by(color, .drop = FALSE) %>% dplyr::summarise(n = dplyr::n(),.groups = "drop") %>% dplyr::filter(color == 'red') %>% dplyr::pull(n)
  message("\t\tThere are ",n_more_percent, "/",nrow(recap)," contrast with more than ",thresholdpct ,"% of genes (significatives and/or DE) lost because of the translation in Entrez ID")
  if(max(recap$pct_sig_not_found, na.rm = TRUE) > max(recap$pct_DE_not_found, na.rm = TRUE)){
    maximum <- max(recap$pct_sig_not_found)
  } else {  
    maximum <- max(recap$pct_DE_not_found)  
  }
  
  g <- ggplot2::ggplot(recap, ggplot2::aes(x = pct_DE_not_found, y = pct_sig_not_found, colour = color)) + 
    ggplot2::geom_point() +
    ggplot2::scale_colour_identity() + 
    ggplot2::facet_grid(control_cluster ~ min.pct + tested_cluster, drop = TRUE) + #en haut: tested_cluster + min.pct; à droite: control_cluster
    ggplot2::annotate("rect", xmin = 0, xmax = thresholdpct, ymin = 0, ymax = thresholdpct, alpha = .2) +
    ggplot2::labs(title="Pourcentages of not found genes", x = "% not found DE genes", y = "% not found signicatives genes") +
    ggplot2::expand_limits(x=c(0,maximum), y=c(0, maximum)) +
    ggplot2::theme(axis.title.x =  ggplot2::element_text(size=9), axis.title.y =  ggplot2::element_text(size=8))
  
  dimension_plot <- dplyr::case_when( length(unique(contrast$tested_cluster)) < 5 ~ 5L,
                                      length(unique(contrast$tested_cluster)) < 50 ~ length(unique(contrast$tested_cluster)),
                                      length(unique(contrast$tested_cluster)) > 50 ~ 49L)
  
  if (DE_type=="1vsAll"){
    ggplot2::ggsave(filename = paste0(DE_type,"_",test,"_pct_not_found.png"), plot = g, width = dimension_plot, height = 1.8, path = output_path_def, )
  }else{
    ggplot2::ggsave(filename = paste0(DE_type,"_",test,"_pct_not_found.png"), plot = g, width = dimension_plot, height = dimension_plot, path = output_path_def)
  }
  
  return(list('recap' = recap,
              'recap_DE_genes' = recap_DE_genes))
}

stat_by_contrast <- function(tested_cluster,control_cluster,min.pct,data, species, pval, thresholdLogFC, thresholdpct, pb = NULL){
  
  Org_db <- get_db(species)['db']
  # Need to use !! to unquote variable name otherwise dplyr::filter will return
  # the all tibble as it will take in priority the name of the columns
  geneList_sig <- data %>% dplyr::filter(tested_cluster == !!tested_cluster & 
                                           control_cluster == !!control_cluster & 
                                           min.pct == !!min.pct) %>% 
    dplyr::filter(adj.P.Val < pval) %>%
    dplyr::select(gene,logFC) %>%
    dplyr::rename(SYMBOL = gene)
  
  #recodage des nom des gènes (ce qui est conseillé)
  #idType("org.Hs.eg.db") #permet de connaitre toutes les conversions possibles
  geneList_sig <- suppressMessages(clusterProfiler::bitr(geneList_sig$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = Org_db)) %>% dplyr::right_join(geneList_sig, by = "SYMBOL")
  
  geneDE <- dplyr::filter(geneList_sig , abs(logFC) > thresholdLogFC)
  
  recap_DE_genes <- tibble::tibble(tested_cluster = tested_cluster, control_cluster = control_cluster, min.pct = min.pct , genes = geneDE$SYMBOL)
  
  All_geneList_sig <- dplyr::filter(geneList_sig,!is.na(ENTREZID)) %>% dplyr::select(ENTREZID,logFC)
  All_geneList <- All_geneList_sig %>% tibble::deframe() %>% sort(decreasing = TRUE)
  All_geneDE <- dplyr::filter(geneDE,!is.na(ENTREZID)) %>% dplyr::select(ENTREZID,logFC) %>% tibble::deframe() %>% sort(decreasing = TRUE)
  
  recap <- tibble::tibble(tested_cluster = tested_cluster,
                          control_cluster = control_cluster,
                          min.pct = min.pct,
                          nb_sig_genes =nrow(geneList_sig),
                          nb_sig_genes_not_found = geneList_sig %>% dplyr::filter(is.na(ENTREZID)) %>% nrow(),
                          nb_sig_genes_used = nrow(All_geneList_sig),
                          nb_DE_genes = nrow(geneDE),
                          nb_DE_genes_not_found = geneDE %>% dplyr::filter(is.na(ENTREZID)) %>% nrow(),
                          nb_DE_genes_used = length(All_geneDE), .rows = 1)
  
  params <- list('tested_cluster' = tested_cluster,
                 'control_cluster' = control_cluster,
                 'min.pct' = min.pct)
  
  if(!is.null(pb)) pb$tick()
  
  return(list("recap" = recap,
              "All_geneList_sig" = All_geneList_sig,
              "All_geneList" = list('gene' = All_geneList,'params' = params),
              "All_geneDE" = list('gene' = All_geneDE, 'params' = params), 
              "recap_DE_genes" = recap_DE_genes))
  
}


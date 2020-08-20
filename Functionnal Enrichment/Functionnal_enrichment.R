suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(purrr))

run_functionnal_enrichment <- function(path, pval = 0.05, thresholdLogFC = 0.58, thresholdpct = 5, species = "mouse", category_msigdb, method = c("WikiPathway", "msigdb", "CellMarkers", "GOterms", "ReactomePA"), gmt_date, gmt_folder, enrichment = TRUE, GSEA = TRUE, comparison = TRUE, verbose = TRUE, force = FALSE) {
  assertthat::assert_that(assertthat::is.dir(path), assertthat::is.writeable(path),assertthat::is.dir(gmt_folder), assertthat::is.readable(gmt_folder))
  assertthat::assert_that(is.numeric(pval), pval >= 0, pval <=1)
  assertthat::assert_that(is.numeric(thresholdpct), thresholdpct >= 0, thresholdpct <=100)
  assertthat::assert_that(is.numeric(thresholdLogFC))
  supported_species <- get_db()
  assertthat::assert_that(is.character(species), is_in(species, supported_species))
  assertthat::assert_that(are_in(method,c("WikiPathway", "msigdb", "CellMarkers", "GOterms", "ReactomePA","DiseaseOntology")))
  assertthat::assert_that(are_in(category_msigdb,c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7")))
  assertthat::assert_that(is.logical(enrichment), is.logical(GSEA), is.logical(comparison), is.logical(verbose), is.logical(verbose))
  
  path <- fs::path_tidy(path)
  
  #Disease ontology warning
  species_syn <- get_db(species)
  if(is.na(species_syn['gmt'])){
    warning("Your species is incompatible with Wiki Pathways")
    method %<>% .[.!="WikiPathway"]
  }
  if(is.na(species_syn['msigdb'])){
    warning("Your species is incompatible with Molecular Signatures (msigdb)")
    method %<>% .[.!="msigdb"]
  }
  if(any(method == "DiseaseOntology") & species != "human"){
    warning("Your species is incompatible with Disease Ontology")
    method %<>% .[.!="DiseaseOntology"]
  }
  
  DE_type <- fs::dir_ls(fs::path(path, "DIFFERENTIAL_EXPRESSION"), type = "directory") %>% fs::path_file()
  output_path <- fs::path(path,"FUNCTIONAL_ENRICHMENT", paste0("pval", pval, "_thresholdLogFC", thresholdLogFC))
  input_path <- fs::path(path,"DE_STAT", paste0("pval", pval, "_thresholdLogFC",thresholdLogFC, "_thresholdpct", thresholdpct))
  # By VS 
  if(verbose) message("Running Functionnal Enrichment : ...")
  purrr::walk(DE_type,run_functionnal_enrichment_by_vs,method = method, species = species, category_msigdb = category_msigdb, gmt_date = gmt_date, gmt_folder = gmt_folder, input_path = input_path, output_path, enrichment = enrichment, GSEA = GSEA, comparison = comparison, verbose = verbose, force = force) 
}

run_functionnal_enrichment_by_vs <- function(DE_type, method, species, category_msigdb, gmt_date, gmt_folder, input_path, output_path, enrichment, GSEA, comparison, verbose, force){
  if(verbose) message("\tDE_type : ",DE_type)
  output_path_Comparisons <- fs::path(output_path, DE_type, "Comparisons")
  test.used <- fs::dir_ls(fs::path(path, "DIFFERENTIAL_EXPRESSION", DE_type), type = "directory", recurse = FALSE) %>% fs::path_file()
  if(verbose) message("\tCreating output folders ...")
  create_folders_enrichment(out.dir = output_path, test = test.used, DE_type = DE_type, enrichment = enrichment, gsea = GSEA)
  
  purrr::walk(test.used, run_functionnal_enrichment_by_test, DE_type = DE_type, method = method, species = species, category_msigdb = category_msigdb, gmt_date = gmt_date, gmt_folder = gmt_folder, input_path = input_path, output_path = output_path, output_path_Comparisons = output_path_Comparisons, enrichment = enrichment, GSEA = GSEA, comparison = comparison, verbose = verbose, force = force)
}

run_functionnal_enrichment_by_test <- function(test, All_geneDE = NULL, All_geneList = NULL, DE_type, method, species, category_msigdb, gmt_date, gmt_folder, input_path, output_path, output_path_Comparisons, enrichment, GSEA, comparison, verbose, force){
  if(verbose) message("\tTest : ", test)
  input_path %<>% fs::path(DE_type, test)
  if (is.null(All_geneDE) & (enrichment | comparison)) All_geneDE <- readRDS(file = fs::path(input_path, paste0(test,"_all_geneDE.rds")))
  if (is.null(All_geneList) & GSEA) All_geneList <- readRDS(file = fs::path(input_path, paste0(test,"_All_geneList.rds")))
  if (any(method == "WikiPathway")){
    gmt_data <- load_wiki_pathway(species = species, gmt_folder, date = gmt_date)
    wpid2gene <- gmt_data[["TERM2GENE"]]
    wpid2name <- gmt_data[["TERM2NAME"]]
  } else {
    wpid2gene <- NULL
    wpid2name <- NULL
  }
  if(any(method == "msigdb")){
    m_t2g <- purrr::set_names(category_msigdb,category_msigdb) %>% purrr::map( ~ load_msigdb(species = species, category = .))
    # m_t2g <- load_msigdb(species = species, category = category_msigdb)
  } else {
    m_t2g <- NULL
  }
  if(any(method == "CellMarkers")){
    cell_markers <- load_cell_markers(species)
  } else {
    cell_markers <- NULL
  }
  
  output_path_def <- fs::path(output_path, DE_type, test)
  
  tested <- test_if_done(test = test, DE_type = DE_type, output_path_def = output_path_def) %>% dplyr::inner_join(tibble::tibble(method = c("enrichment", "GSEA", "comparisons"), user = c(enrichment, GSEA, comparison)), by = "method")
  if(verbose){
    pb_enrich <- progress::progress_bar$new(format = "\t\tComputing : |:bar| :percent (:elapsedfull)",total = length(All_geneDE), complete = "█", incomplete = " ", current = " ", clear = FALSE)
    pb_gsea <- progress::progress_bar$new(format = "\t\tComputing : |:bar| :percent (:elapsedfull)",total = length(All_geneList), complete = "█", incomplete = " ", current = " ", clear = FALSE)
    pb_plot_enrich <- progress::progress_bar$new(format = "\t\tPlotting : |:bar| :percent (:elapsedfull)",total = length(All_geneDE), complete = "█", incomplete = " ", current = " ", clear = FALSE)
    pb_plot_gsea <- progress::progress_bar$new(format = "\t\tPlotting : |:bar| :percent (:elapsedfull)",total = length(All_geneList), complete = "█", incomplete = " ", current = " ", clear = FALSE)
  } else {
    pb_enrich <- NULL
    pb_gsea <- NULL
    pb_plot_enrich <- NULL
    pb_plot_gsea <- NULL
  }

  temp <- tested %>% 
    dplyr::filter(method == "enrichment") 
  if(force | (temp$user & !temp$data & !temp$Plot)){
    if(verbose) {
      message("\t\tEnrichment Analysis ...")
      pb_enrich$tick(0)
    }
    purrr::map(All_geneDE, run_enrichment, method = method, category_msigdb = category_msigdb, wpid2gene = wpid2gene, wpid2name = wpid2name, m_t2g = m_t2g, cell_markers = cell_markers, species = species, output_path_def = output_path_def, DE_type = DE_type, pb = pb_enrich) %>%
      result_message_summary(path = fs::path(output_path_def,"Enrichment_plot")) %T>% 
      saveRDS(file = fs::path(output_path_def, paste0(test,"_enrichment.rds"))) %>%
      purrr::walk(plot_all_enrichment, output_path_def = output_path_def, GSEA = FALSE,pb = pb_plot_enrich)
  } else if(temp$user & temp$data & !temp$Plot){
    message("\t\tEnrichment : Data already exists but plot are incomplete")
    # if(verbose) pb_plot_enrich$tick(0)
    readRDS(fs::path(output_path_def, paste0(test,"_enrichment.rds"))) %>% 
      purrr::walk(plot_all_enrichment, output_path_def = output_path_def, GSEA = FALSE, pb = pb_plot_enrich)
  } else if(temp$user & temp$data & temp$Plot){
    message("\t\tEnrichment : Data and Plots already exist, skipping ... ")
  }
  gc(verbose = FALSE)
  
  temp <- tested %>% 
    dplyr::filter(method == "GSEA")
  if(force | (temp$user & !temp$data & !temp$Plot)){
    if(verbose) {
      message("\t\tGSEA Analysis ...")
      pb_gsea$tick(0) 
    }
    purrr::map(All_geneList, run_GSEA, method = method, category_msigdb = category_msigdb, wpid2gene = wpid2gene, wpid2name = wpid2name, m_t2g = m_t2g, cell_markers = cell_markers, species = species, output_path_def = output_path_def, DE_type = DE_type, pb = pb_gsea) %>%
      result_message_summary(path = fs::path(output_path_def,"GSEA_plot")) %T>% 
      saveRDS(file = fs::path(output_path_def, paste0(test,"_gsea.rds"))) %>% 
      purrr::walk(plot_all_enrichment, output_path_def = output_path_def, GSEA = TRUE,pb = pb_plot_gsea)
  } else if(temp$user & temp$data & !temp$Plot){
    message("\t\tGSEA : Data already exists but plot are incomplete")
    if(verbose) pb_plot_gsea$tick(0)
    readRDS(fs::path(output_path_def, paste0(test,"_gsea.rds"))) %>% 
      purrr::walk(plot_all_enrichment, output_path_def = output_path_def, GSEA = TRUE, pb = pb_plot_gsea)
  } else if(temp$user & temp$data & temp$Plot){
    message("\t\tGSEA : Data and Plots already exist, skipping ... ")
  }
  
  rm(All_geneList,wpid2gene,wpid2name,m_t2g,cell_markers)
  gc(verbose = FALSE)
  
  temp <- tested %>% 
    dplyr::filter(method == "comparisons")

  if(length(All_geneDE) > 2 & comparison & any(!temp$data,!temp$Plot)){
    if(verbose) message("\t\tComparison Analysis...")
    cluster_comparison(All_geneDE = All_geneDE, output_path_Comparisons = output_path_Comparisons, test = test, DE_type = DE_type, fun = "enrichGO", species = species)
  } else {
    message("\t\tComparison Analysis already performed, skipping ...")
  }
}

run_enrichment <- function(input, method, category_msigdb = NULL, wpid2gene = NULL, wpid2name = NULL, m_t2g = NULL, cell_markers = NULL, species, output_path_def, DE_type, pb = NULL){
  Org_db = get_db(species)['db']
  if( any(method == "WikiPathway") & (is.null(wpid2name) | is.null(wpid2gene)) ){
    stop("To run WikiPathway please specify `wpid2gene`, `wpid2name` and `Org_db` ")
  }
  if(any(method == "msigdb") & is.null(m_t2g)){
    stop("To run msigdb please specify `m_t2g` and `Org_db`")
  }
  if(any(method == "CellMarkers") & is.null(cell_markers)){
    stop("To run CellMarkers please specify `cell_markers` and `Org_db`")
  }
  if(any(method == "DiseaseOntology") & Org_db !=  "org.Hs.eg.db"){
    warning("Disease ontology only works for human (org.Hs.eg.db) ")
    method <- method[!which(method == DiseaseOntology)]
  }
  
  out <- list()
  out[['result']] <- list()
  out[['messages']] <- list()
  out[['warnings']] <- list()
  
  if(any(method == "WikiPathway")) {
    enrich_wp <-  purrr::quietly(WikiPathway_enrichment)(genes = names(input$gene), wpid2gene, wpid2name, Org_db, tested_cluster = input$params$tested_cluster, control_cluster = input$params$control_cluster, min.pct = input$params$min.pct, output_path_def, DE_type )
    out[['result']][['Wiki Pathway']] <- purrr::pluck(enrich_wp, "result")
    out[['warnings']][['Wiki Pathway']] <- purrr::pluck(enrich_wp, "warnings")
    out[['messages']][['Wiki Pathway']] <- purrr::pluck(enrich_wp, "messages")
  }
  if(any(method == "msigdb")) {
    enrich_msig <-  purrr::map(m_t2g, ~ purrr::quietly(msigdb_enrichment)(genes = names(input$gene),m_t2g = ., Org_db = Org_db, tested_cluster = input$params$tested_cluster, control_cluster = input$params$control_cluster, min.pct = input$params$min.pct, output_path_def = output_path_def, DE_type = DE_type))
    enrich_msig %<>%
      purrr::set_names(paste0('Molecular Signatures (', names(.), ')')) 
    out <- assing_out_values(enrich_msig,out,"result")
    out <- assing_out_values(enrich_msig,out,"messages")
    out <- assing_out_values(enrich_msig,out,"warnings")
  }
  if(any(method == "CellMarkers")) {
    enrich_cm <-  purrr::quietly(CellMarkers_enrichment)(names(input$gene),cell_markers, Org_db, tested_cluster = input$params$tested_cluster, control_cluster = input$params$control_cluster, min.pct = input$params$min.pct, output_path_def, DE_type)
    out[['result']][['Cell Markers']] <- purrr::pluck(enrich_cm, "result")
    out[['warnings']][['Cell Markers']] <- purrr::pluck(enrich_cm, "warnings")
    out[['messages']][['Cell Markers']] <- purrr::pluck(enrich_cm, "messages")
  }
  if(any(method == "DiseaseOntology")) {
    enrich_do <-  purrr::quietly(DiseaseOntology_enrichment)(names(input$gene), tested_cluster = input$params$tested_cluster, control_cluster = input$params$control_cluster, min.pct = input$params$min.pct, output_path_def, DE_type)
    out[['result']][['Disease Ontology']] <- purrr::pluck(enrich_do, "result")
    out[['warnings']][['Disease Ontology']] <- purrr::pluck(enrich_do, "warnings")
    out[['messages']][['Disease Ontology']] <- purrr::pluck(enrich_do, "messages")
  }
  if(any(method == "GOterms")) {
    enrich_go <-  purrr::quietly(GOterms_enrichment)(names(input$gene), Org_db, tested_cluster = input$params$tested_cluster, control_cluster = input$params$control_cluster, min.pct = input$params$min.pct, output_path_def, DE_type)
    out[['result']][['Gene Ontology']] <- purrr::pluck(enrich_go, "result")
    out[['warnings']][['Gene Ontology']] <- purrr::pluck(enrich_go, "warnings")
    out[['messages']][['Gene Ontology']] <- purrr::pluck(enrich_go, "messages")
  }
  if(any(method == "ReactomePA")){
    enrich_rPA <-  purrr::quietly(Reactome_enrichment)(names(input$gene), species, tested_cluster = input$params$tested_cluster, control_cluster = input$params$control_cluster, min.pct = input$params$min.pct, output_path_def, DE_type)
    out[['result']][['Reactome Pathways Analysis']] <- purrr::pluck(enrich_rPA, "result")
    out[['warnings']][['Reactome Pathways Analysis']] <- purrr::pluck(enrich_rPA, "warnings")
    out[['messages']][['Reactome Pathways Analysis']] <- purrr::pluck(enrich_rPA, "messages")
  }
  
  out[['result']][['params']] <- list('tested_cluster' = input$params$tested_cluster,
                          'control_cluster' = input$params$control_cluster,
                          'min.pct' = input$params$min.pct)
  out[['result']][['gene']] <- input$gene
  if(!is.null(pb)) pb$tick()
  return(out)
}

run_GSEA <- function(input, method, category_msigdb = NULL, wpid2gene = NULL, wpid2name = NULL, m_t2g = NULL, cell_markers = NULL, species, output_path_def, DE_type, pb = NULL){
  Org_db = get_db(species)['db']
  if( any(method == "WikiPathway") & (is.null(wpid2name) | is.null(wpid2gene)) ){
    stop("To run WikiPathway please specify `wpid2gene`, `wpid2name` and `Org_db` ")
  }
  if(any(method == "msigdb") & is.null(m_t2g) ){
    stop("To run msigdb please specify `m_t2g` and `Org_db`")
  }
  if(any(method == "CellMarkers") & is.null(cell_markers)){
    stop("To run CellMarkers please specify `cell_markers` and `Org_db`")
  }
  if(any(method == "DiseaseOntology") & Org_db !=  "org.Hs.eg.db"){
    warning("Disease ontology only works for human (org.Mm.eg.db) ")
    method <- method[!which(method == DiseaseOntology)]
  }
  
  out <- list()
  out[['result']] <- list()
  out[['messages']] <- list()
  out[['warnings']] <- list()

  
  if(any(method == "WikiPathway")) {
    gsea_wp <- purrr::quietly(WikiPathway_GSEA)(genes = input$gene,wpid2gene,wpid2name,Org_db, tested_cluster = input$params$tested_cluster, control_cluster = input$params$control_cluster, min.pct = input$params$min.pct, output_path_def, DE_type )
    out[['result']][['Wiki Pathway']] <- purrr::pluck(gsea_wp, "result")
    out[['warnings']][['Wiki Pathway']] <- purrr::pluck(gsea_wp, "warnings")
    out[['messages']][['Wiki Pathway']] <- purrr::pluck(gsea_wp, "messages")
  }
  if(any(method == "msigdb")) {
    gsea_msig <-  purrr::map(m_t2g, ~ purrr::quietly(msigdb_GSEA)(genes = input$gene, m_t2g = ., Org_db = Org_db, tested_cluster = input$params$tested_cluster, control_cluster = input$params$control_cluster, min.pct = input$params$min.pct, output_path_def = output_path_def, DE_type = DE_type))
    gsea_msig %<>%
      purrr::set_names(paste0('Molecular Signatures (', names(.), ')')) 
    out <- assing_out_values(gsea_msig,out,"result")
    out <- assing_out_values(gsea_msig,out,"messages")
    out <- assing_out_values(gsea_msig,out,"warnings")
  }
  if(any(method == "DiseaseOntology")) {
    gsea_do <- purrr::quietly(DiseaseOntology_GSEA)(genes = input$gene, Org_db,tested_cluster = input$params$tested_cluster, control_cluster = input$params$control_cluster, min.pct = input$params$min.pct, output_path_def, DE_type)
    out[['result']][['Disease Ontology']] <- purrr::pluck(gsea_do, "result")
    out[['warnings']][['Disease Ontology']] <- purrr::pluck(gsea_do, "warnings")
    out[['messages']][['Disease Ontology']] <- purrr::pluck(gsea_do, "messages")
  }
  if(any(method == "GOterms")) {
    gsea_go <- purrr::quietly(GOterm_GSEA)(genes = input$gene, Org_db,tested_cluster = input$params$tested_cluster, control_cluster = input$params$control_cluster, min.pct = input$params$min.pct, output_path_def, DE_type)
    out[['result']][['Gene Ontology']] <- purrr::pluck(gsea_go, "result")
    out[['warnings']][['Gene Ontology']] <- purrr::pluck(gsea_go, "warnings")
    out[['messages']][['Gene Ontology']] <- purrr::pluck(gsea_go, "messages")
  }
  if(any(method == "ReactomePA")){
    gsea_rPA <- purrr::quietly(Reactome_GSEA)(input$gene, species, tested_cluster = input$params$tested_cluster, control_cluster = input$params$control_cluster, min.pct = input$params$min.pct, output_path_def, DE_type)
    out[['result']][['Reactome Pathways Analysis']] <- purrr::pluck(gsea_rPA, "result")
    out[['warnings']][['Reactome Pathways Analysis']] <- purrr::pluck(gsea_rPA, "warnings")
    out[['messages']][['Reactome Pathways Analysis']] <- purrr::pluck(gsea_rPA, "messages")
  }
  out[['result']][['params']] <- list('tested_cluster' = input$params$tested_cluster,
                          'control_cluster' = input$params$control_cluster,
                          'min.pct' = input$params$min.pct)
  out[['result']][['gene']] <- input$gene
  if(!is.null(pb)) pb$tick()
  return(out)
}

WikiPathway_enrichment <- function(genes,wpid2gene,wpid2name,Org_db, tested_cluster, control_cluster, min.pct, output_path_def, DE_type ){
  enrich <- clusterProfiler::enricher(genes, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)

  if(!is.null(enrich) && NROW(enrich) > 0){
    enrich <- clusterProfiler::setReadable(enrich, OrgDb = Org_db, keyType = "ENTREZID") #mettre les noms de gènes à la place de l'ID
    enrich %>%
      slot("result") %>% # extract result slot
      tibble::as_tibble() %>% 
      tibble::add_column(tested_cluster = tested_cluster, control_cluster = control_cluster, min.pct = min.pct, .before = 1) %>%
      write_enrichment_results("WikiPathways_enrichment", output_path_def, DE_type)
    return(enrich)
  } else {
    return(NULL)
  }
}

WikiPathway_GSEA <- function(genes,wpid2gene,wpid2name,Org_db, tested_cluster, control_cluster, min.pct, output_path_def, DE_type ){
  gsea <- purrr::possibly(clusterProfiler::GSEA, otherwise = NULL)(genes, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=FALSE)
  if(!is.null(gsea) && NROW(gsea) > 0){
    gsea <- clusterProfiler::setReadable(gsea, OrgDb = Org_db, keyType = "ENTREZID") #mettre les noms de gènes à la place de l'ID
    gsea %>% 
      slot("result") %>% # extract result slot
      tibble::as_tibble() %>% 
      tibble::add_column(tested_cluster = tested_cluster, control_cluster = control_cluster, min.pct = min.pct, .before = 1) %>%
      write_enrichment_results("WikiPathways_GSEA",output_path_def, DE_type)
    return(gsea)
  }  else {
    return(NULL)
  }
  
}

msigdb_enrichment <- function(genes,m_t2g,Org_db, tested_cluster, control_cluster, min.pct, output_path_def, DE_type){
  enrich <- clusterProfiler::enricher(genes, TERM2GENE=m_t2g)
  if(!is.null(enrich) && NROW(enrich) > 0){
    enrich <- clusterProfiler::setReadable(enrich, OrgDb = Org_db, keyType = "ENTREZID") #mettre les noms de gènes à la place de l'ID
    enrich %>%
      slot("result") %>% # extract result slot
      tibble::as_tibble() %>% 
      tibble::add_column(tested_cluster = tested_cluster, control_cluster = control_cluster, min.pct = min.pct, .before = 1) %>%
      write_enrichment_results("msigdb_enrichment",output_path_def, DE_type )
    return(enrich)
  } else {
    return(NULL)
  }
}

msigdb_GSEA <- function(genes,m_t2g,Org_db, tested_cluster, control_cluster, min.pct, output_path_def, DE_type){
  gsea <- purrr::possibly(clusterProfiler::GSEA,otherwise = NULL)(genes, TERM2GENE=m_t2g, verbose = FALSE)
  if(!is.null(gsea) && NROW(gsea) > 0){
    gsea <- clusterProfiler::setReadable(gsea, OrgDb = Org_db, keyType = "ENTREZID") #mettre les noms de gènes à la place de l'ID
      
    gsea %>%
      slot("result") %>% # extract result slot
      tibble::as_tibble() %>% 
      tibble::add_column(tested_cluster = tested_cluster, control_cluster = control_cluster, min.pct = min.pct, .before = 1) %>%
      write_enrichment_results("msigdb_GSEA" ,output_path_def, DE_type)
    return(gsea)
  } else {
    return(NULL)
  }
}

CellMarkers_enrichment <- function(genes,cell_markers,Org_db, tested_cluster, control_cluster, min.pct, output_path_def, DE_type){
  enrich <- clusterProfiler::enricher(genes, TERM2GENE=cell_markers)  
  if(!is.null(enrich) && NROW(enrich) > 0){
    enrich <- clusterProfiler::setReadable(enrich, OrgDb = Org_db, keyType = "ENTREZID") #mettre les noms de gènes à la place de l'ID
    enrich %>%
      slot("result") %>% # extract result slot
      tibble::as_tibble() %>% 
      tibble::add_column(tested_cluster = tested_cluster, control_cluster = control_cluster, min.pct = min.pct, .before = 1) %>%
      write_enrichment_results("CellMarkers_enrichment",output_path_def, DE_type )
    return(enrich)
  } else {
    return(NULL)
  }
}

DiseaseOntology_enrichment <- function(genes, tested_cluster, control_cluster, min.pct, output_path_def, DE_type){
  enrich <- DOSE::enrichDO(gene          = genes,
                           ont           = "DO",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH",
                           universe      = names(genes),
                           minGSSize     = 5,
                           maxGSSize     = 500,
                           qvalueCutoff  = 0.05,
                           readable      = TRUE)
  if(!is.null(enrich) && NROW(enrich) > 0){
    tibble::as_tibble(enrich@result) %>% 
      tibble::add_column(tested_cluster = tested_cluster, control_cluster = control_cluster, min.pct = min.pct, .before = 1) %>%
      write_enrichment_results("DiseaseOntology_enrichment",output_path_def, DE_type)
    return(enrich)
  } else {
    return(NULL)
  }
}

DiseaseOntology_GSEA <- function(genes, Org_db,tested_cluster, control_cluster, min.pct, output_path_def, DE_type){
  gsea <- DOSE::gseDO(geneList      = genes,
                      nPerm         = 100,
                      minGSSize     = 120,
                      pvalueCutoff  = 0.2,
                      pAdjustMethod = "BH",
                      verbose       = FALSE) 
  if(!is.null(gsea) && NROW(gsea) > 0){
    gsea <- clusterProfiler::setReadable(gsea, OrgDb = Org_db, keyType = "ENTREZID") #mettre les noms de gènes à la place de l'ID
    
    gsea %>% 
      slot("result") %>% # extract result slot
      tibble::as_tibble() %>% 
      tibble::add_column(tested_cluster = tested_cluster, control_cluster = control_cluster, min.pct = min.pct, .before = 1) %>%
      write_enrichment_results("DiseaseOntology_GSEA",output_path_def, DE_type )
    return(gsea)
  } else {
    return(NULL)
  }
}

GOterms_enrichment <- function(genes, Org_db, tested_cluster, control_cluster, min.pct, output_path_def, DE_type){
  enrich <- clusterProfiler::enrichGO(gene          = genes,
                                      OrgDb         = Org_db,
                                      ont           = "BP",
                                      pvalueCutoff  = 0.01,
                                      pAdjustMethod = "BH",
                                      qvalueCutoff  = 0.05,
                                      readable      = TRUE)
  if(!is.null(enrich) && NROW(enrich) > 0){
    tibble::as_tibble(enrich@result) %>% 
      tibble::add_column(tested_cluster = tested_cluster, control_cluster = control_cluster, min.pct = min.pct, .before = 1) %>% 
      write_enrichment_results("GOterms_enrichment",output_path_def, DE_type )
    return(enrich)
  } else {
    return(NULL)
  }
}

GOterm_GSEA <- function(genes, Org_db,tested_cluster, control_cluster, min.pct, output_path_def, DE_type){
  gsea <- purrr::possibly(clusterProfiler::gseGO,otherwise = NULL)(geneList      = genes,
                                 OrgDb         = Org_db,
                                 ont           = "BP",
                                 nPerm         = 1000,
                                 minGSSize     = 120,
                                 pvalueCutoff  = 0.05,
                                 pAdjustMethod = "BH",
                                 verbose       = FALSE) 
  if(!is.null(gsea) && NROW(gsea) > 0){
   gsea <- clusterProfiler::setReadable(gsea, OrgDb = Org_db, keyType = "ENTREZID") #mettre les noms de gènes à la place de l'ID
   
   gsea %>%
      slot("result") %>% # extract result slot
      tibble::as_tibble() %>% 
      tibble::add_column(tested_cluster = tested_cluster, control_cluster = control_cluster, min.pct = min.pct, .before = 1) %>%
      write_enrichment_results("GOterms_GSEA",output_path_def, DE_type )
    return(gsea)
  } else {
    return(NULL)
  }
}

Reactome_enrichment <- function(genes, species, tested_cluster, control_cluster, min.pct, output_path_def, DE_type){
  enrich <- ReactomePA::enrichPathway(gene = genes, organism = species, pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, readable = TRUE) 
  if(!is.null(enrich) && NROW(enrich) > 0){
    tibble::as_tibble(enrich@result) %>% 
      tibble::add_column(tested_cluster = tested_cluster, control_cluster = control_cluster, min.pct = min.pct, .before = 1) %>% 
      write_enrichment_results("Reactome_enrichment",output_path_def, DE_type )
    return(enrich)
  } else {
    return(NULL)
  }
}

Reactome_GSEA <- function(genes, species, tested_cluster, control_cluster, min.pct, output_path_def, DE_type){
  Org_db <- get_db(organism = species)['db']
  gsea <- purrr::possibly(ReactomePA::gsePathway, otherwise = NULL)(geneList = genes, organism = species, pvalueCutoff = 0.05, pAdjustMethod = "BH", verbose = FALSE)
  if(!is.null(gsea) && NROW(gsea) > 0){
    gsea <- clusterProfiler::setReadable(gsea, OrgDb = Org_db, keyType = "ENTREZID")#mettre les noms de gènes à la place de l'ID
    gsea %>%
      slot("result") %>% # extract result slot
      tibble::as_tibble() %>% 
      tibble::add_column(tested_cluster = tested_cluster, control_cluster = control_cluster, min.pct = min.pct, .before = 1) %>% 
      write_enrichment_results("Reactome_GSEA" ,output_path_def, DE_type)
    return(gsea)
  } else {
    return(NULL)
  }
}

write_enrichment_results <- function(df,filename, output_path_def, DE_type){
  if(fs::file_exists(fs::path(output_path_def,paste0(filename,"_",DE_type,".txt")))){
    readr::write_tsv(df, path = fs::path(output_path_def,paste0(filename,"_",DE_type,".txt")), append = TRUE, col_names = FALSE, na = "NA")
  } else {
    readr::write_tsv(df, path = fs::path(output_path_def,paste0(filename,"_",DE_type,".txt")), col_names = TRUE, append = FALSE, na = "NA")
  }
}

create_folders_enrichment <- function(out.dir, test, DE_type, enrichment = TRUE, gsea = TRUE){
  
  if(!enrichment & !gsea){
    stop("No directory to create")
  }
  
  dir <- c("Enrichment_plot/dotplot","Enrichment_plot/upsetplot","Enrichment_plot/emapplot","Enrichment_plot/goplot",
           "Enrichment_plot/Pathways","GSEA_plot/dotplot","GSEA_plot/emapplot","GSEA_plot/gseaplot",
           "GSEA_plot/upsetplot","GSEA_plot/Pathways")
  
  paths <- tidyr::expand_grid(out.dir,DE_type,test,dir) %>% 
    dplyr::transmute(path = stringr::str_c(out.dir,DE_type,test,dir, sep = "/")) %>% 
    dplyr::pull(path) 
  
  if(!enrichment) paths %<>% stringr::str_subset("GSEA")
  if(!gsea) paths %<>% stringr::str_subset("Enrichment")
  
  tidyr::expand_grid(out.dir, DE_type, dir = "Comparisons") %>% 
    dplyr::transmute(path = stringr::str_c(out.dir,DE_type,dir, sep = "/")) %>% 
    dplyr::pull(path) %>%
    vctrs::vec_c(paths) %>% 
    fs::dir_exists() %>% 
    .[. == FALSE] %>% 
    names() %>%
    fs::dir_create(recurse = TRUE)
}

plot_all_enrichment <- function(input, output_path_def, GSEA, pb = NULL){
  params <- input$params
  input$params <- NULL
  logFC <- input$gene
  input$gene <- NULL
  enrich_method <- ifelse(GSEA,"gsea","enrich")
  purrr::iwalk(input, plot_enrichment, params = params, output_path_def = output_path_def, enrich_method = enrich_method, logFC = logFC)
  if(!is.null(pb)) pb$tick()
}

plot_enrichment <- function(enrich, method, params, enrich_method, output_path_def, logFC){
  output_path_plot <- switch(enrich_method,
                             gsea = fs::path(output_path_def,'GSEA_plot'),
                             enrich = fs::path(output_path_def,'Enrichment_plot'))

  if(class(enrich) == "gseaResult"){
    gene_logFC <- enrich@gene2Symbol %>% 
      tibble::enframe(name = "ENTREZID", value = "gene") %>% 
      dplyr::inner_join(tibble::enframe(enrich@geneList, name = "ENTREZID", value = "logFC"), by = "ENTREZID")
  }
  
  p1 <- clusterProfiler::dotplot(enrich, showCategory=20) +
    ggplot2::labs(title=paste0("Dotplot for ",method,": cluster ", params$tested_cluster, " vs ", params$control_cluster, " for min.pct ", params$min.pct)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=8), axis.text.y = ggplot2::element_text(size=10), axis.title.x = ggplot2::element_text(size = 10), plot.title.position = 'plot')
  
  utils <- p1$data %>% 
    dplyr::select(Description,x) %>% 
    dplyr::mutate(y_len = stringr::str_length(Description)) %>% 
    dplyr::summarise(max_y = max(y_len), min_x = min(x), max_x = max(x), .groups = "drop")
  
  n_breaks <- ggplot2::ggplot_build(p1)$layout$panel_params[[1]]$x$breaks %>% length() - 1
  # n_breaks <- labeling::extended(range(p$data$x)[1], range(p$data$x)[2], m = 5 ) %>% length() - 1 
  width = utils$max_y * 0.07 + n_breaks * 1.15 + 1.25
  suppressMessages(ggplot2::ggsave(filename = paste0(method,": cluster ", params$tested_cluster, " vs ", params$control_cluster, " for min.pct ", params$min.pct,"_dotplot.png"), plot = p1, path = fs::path(output_path_plot,"dotplot"), width = width, height = width / 1.3))
  
  # UPSET PLOT
  p_temp <- enrichplot::upsetplot(enrich, n=10) 
  # Workaround to only display 20 first intersection in upset plot
  # as n_intersections in scale_x_upset gives an error.
  
  p_temp$data %>% 
    dplyr::select(Description) %>%
    dplyr::mutate(Description = as.character(Description)) %>% 
    dplyr::group_by(Description) %>% 
    dplyr::add_count(sort = TRUE) %>% 
    dplyr::distinct() %>% 
    dplyr::ungroup() %>% 
    dplyr::slice(1:20) %>% 
    dplyr::pull(Description) -> temp
  
  p2 <- p_temp$data %>% 
    dplyr::filter(Description %in% temp) %>% 
    purrr::when(enrich_method == "gsea" ~ dplyr::inner_join(.,gene_logFC, by = "gene") %>% ggplot2::ggplot(ggplot2::aes(x = Description, y = logFC)) + ggplot2::geom_boxplot() + ggplot2::geom_jitter(width = .2, alpha = .6) + ggplot2::ylab("logFC"),
                enrich_method == "enrich" ~ ggplot2::ggplot(.,ggplot2::aes(x = Description)) + ggplot2::geom_bar() + ggplot2::ylab("Number of genes in common"))
  p2 <- p2 + 
    ggupset::scale_x_upset() + ggplot2::xlab(NULL) + 
    ggplot2::labs(title=paste0("Upsetplot for ",method,": cluster ", params$tested_cluster, " vs ", params$control_cluster, " for min.pct ", params$min.pct)) +
    ggplot2::theme(plot.title = ggplot2::element_text(size=10), axis.title.y = ggplot2::element_text(size=10))
  
  if(enrich_method == 'enrich'){
    p2 <- p2 + ggplot2::geom_text(stat='count', ggplot2::aes(label=..count..), vjust = -0.5) + 
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .15)))
    if(method == "Gene Ontology"){
      suppressMessages(p <- enrichplot::goplot(enrich, showCategory = 5) +
        ggplot2::labs(title = paste0("GO plot for ",method,": cluster ", params$tested_cluster, " vs ", params$control_cluster, " for min.pct ", params$min.pct)))
      suppressMessages(ggplot2::ggsave(filename = paste0(method,": cluster ", params$tested_cluster, " vs ", params$control_cluster, " for min.pct ", params$min.pct,"_goplot.png"), plot = p, width=8, height=5, path = fs::path(output_path_plot,"goplot")))
    }
  } else {
    p4 <- enrichplot::gseaplot2(enrich, geneSetID = 1: dplyr::case_when(length(enrich$Description) < 10 ~ length(enrich$Description),
                                                                        TRUE ~ 10L)) +
      ggplot2::labs(title = paste0("GSEA plot for ",method,": cluster ", params$tested_cluster, " vs ", params$control_cluster, " for min.pct ", params$min.pct))
    suppressMessages(ggplot2::ggsave(filename = paste0(method,": cluster ", params$tested_cluster, " vs ", params$control_cluster, " for min.pct ", params$min.pct,"_gseaplot.png"), plot = p4, width=10, height=11, path = fs::path(output_path_plot,"gseaplot")))
  }
  
  suppressMessages(ggplot2::ggsave(filename = paste0(method,": cluster ", params$tested_cluster, " vs ", params$control_cluster, " for min.pct ", params$min.pct,"_upsetplot.png"), plot = p2, width=10, height=6, path = fs::path(output_path_plot,"upsetplot")))
  
  suppressMessages(p3 <- enrichplot::emapplot(enrich, showCategory=20) +
    ggplot2::labs(title = paste0("Enrichment Map for ",method,": cluster ", params$tested_cluster, " vs ", params$control_cluster, " for min.pct ", params$min.pct)))
  suppressMessages(ggplot2::ggsave(filename = paste0(method,": cluster ", params$tested_cluster, " vs ", params$control_cluster, " for min.pct ", params$min.pct,"_emapplot.png"), plot = p3, width=10, height=10, path = fs::path(output_path_plot,"emapplot")))
  
  if(method == "Reactome Pathways Analysis"){
    enrich@result %>% 
      dplyr::slice(1:10) %>% 
      dplyr::pull(Description) %>% 
      purrr::set_names(.) %>% 
      purrr::map(purrr::possibly(viewPathway_igraph, otherwise = NULL), organism = "mouse", readable = TRUE, foldChange = logFC) %>% 
      purrr::compact() -> igraph_pathways 
    
    purrr::imap(igraph_pathways, plot_igraph_static, params = params) %>%  
      saveRDS(file = fs::path(output_path_plot,"Pathways",paste0("cluster ", params$tested_cluster, " vs ", params$control_cluster, " for min.pct ", params$min.pct,"_dotplot.rds")))
    
    comp <- paste0(params$tested_cluster, " vs ", params$control_cluster, " with min.pct ", params$min.pct)
    fs::path(output_path_plot,"Pathways", comp) %>% 
      fs::dir_exists() %>% 
      .[. == FALSE] %>% 
      names() %>%
      fs::dir_create()
    # https://github.com/ramnathv/htmlwidgets/issues/299
    # Issue in saving if not in the current directory
    # Use of possbly beacuse of an error with pathway : "Carboxyterminal post-translational modifications of tubulin"
    # invalid vertex id 
    purrr::imap(igraph_pathways, purrr::possibly(plot_igraph_interactive, otherwise = NULL), params = params) %>%
      purrr::compact() %>%
      purrr::iwalk(~ withr::with_dir(fs::path(output_path_plot,"Pathways",comp), visNetwork::visSave(.x, file = paste0(stringr::str_replace_all(.y,'/','-'),".html"))))
  }
}

cluster_comparison <- function(All_geneDE, output_path_Comparisons, test, DE_type, fun = "enrichGO",species){
  Org_db <- get_db(species)['db']
  data <- purrr::map_df(All_geneDE, ~ tibble::enframe(.x[['gene']]) %>% 
                          tibble::add_column(group = paste0(.x[['params']][["tested_cluster"]],"vs", .x[['params']][["control_cluster"]]), 
                                             min.pct = factor(.x[['params']][["min.pct"]] * 100, levels = sort(unique(.x[['params']][["min.pct"]] * 100)))),.name_repair = 'minimal') %>% 
    dplyr::rename(ENTREZID = name, logFC = value) %>%
    dplyr::select(ENTREZID,group, min.pct)
  cluster_comp <- purrr::possibly(clusterProfiler::compareCluster, otherwise = NULL)(ENTREZID ~ group + min.pct, fun = "enrichGO", OrgDb = Org_db, data = data)
  if(!is.null(cluster_comp)){
    saveRDS(cluster_comp, file = fs::path(output_path_Comparisons, paste0(test,"_",DE_type,"_cluster_comaprison.rds")))
    if(length(unique(dplyr::pull(data,group))) > 20){
      unique(stringr::str_extract(cluster_comp@compareClusterResult$group,"^[0-9]+vs")) %>% 
        stringr::str_c("^",.) %>% 
        purrr::walk(plot_many_comp, cluster_comp = cluster_comp, output_path_Comparisons = output_path_Comparisons, test = test, DE_type = DE_type)
    } else {
      cluster_comp@compareClusterResult %<>% dplyr::mutate(group = factor(group, stringr::str_sort(unique(group), numeric = TRUE)))
      p <- enrichplot::dotplot(cluster_comp,x=~group, showCategory = 5) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size=8, angle = 90), axis.text.y = ggplot2::element_text(size=8), axis.title.x = ggplot2::element_text(size = 8)) + 
        ggplot2::facet_grid(~min.pct, labeller = ggplot2::labeller(min.pct = label_pct)) + ggplot2::ggtitle("Comparison for all cluster ")
      ggplot2::ggsave(filename = paste0(test,"_",DE_type,"_cluster_all_dotplot.png"), plot = p, width=15, height=6, path = output_path_Comparisons)
    }
  }
}

label_pct <- function(value){
  out <- stringr::str_c("Percentage min : ",value," %")
}

plot_many_comp <- function(comp,cluster_comp,output_path_Comparisons,test,DE_type){
  cluster_comp@compareClusterResult %<>% dplyr::filter(stringr::str_detect(group,comp)) %>% dplyr::mutate(group = factor(group, stringr::str_sort(unique(group), numeric = TRUE)))
  p <- enrichplot::dotplot(cluster_comp,x=~group, showCategory = 3) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=8, angle = 90), axis.text.y = ggplot2::element_text(size=8), axis.title.x = ggplot2::element_text(size = 8)) + ggplot2::facet_grid(~min.pct, labeller = ggplot2::labeller(min.pct = label_pct)) + ggplot2::ggtitle(paste0("Comparison for cluster ",stringr::str_extract(comp, "[0-9]+")))
  ggplot2::ggsave(filename = paste0(test,"_",DE_type,"_cluster_",stringr::str_extract(comp, "[0-9]+"),"_dotplot.png"), plot = p, width=12, height=5, path = output_path_Comparisons)
}

get_db <- function(organism) {
  supported_species <- c("anopheles","arabidopsis","bovine","canine","celegans","chicken","chimp","coelicolor","ecolik12","ecsakai","fly","gondii","human","malaria","mouse","pig","rat","rhesus","xenopus","yeast","zebrafish"
)
  if (rlang::is_missing(organism)){
    # return list of supported organism if no species specified
    supported_species
  } else {
    assertthat::assert_that(is_in(organism, supported_species))
    annoDb <- switch(organism,
                     anopheles   = c( db = "org.Ag.eg.db", gmt = "Anopheles_gambiae", msigdb = NA),
                     arabidopsis = c( db = "org.At.tair.db", gmt = "Arabidopsis_thaliana", msigdb = NA),
                     bovine      = c( db = "org.Bt.eg.db", gmt = "Bos_taurus", msigdb = "Bos_taurus"),
                     canine      = c( db = "org.Cf.eg.db", gmt = "Canis_familiaris", msigdb = "Canis lupus familiaris"),
                     celegans    = c( db = "org.Ce.eg.db", gmt = "Caenorhabditis_elegans", msigdb = "Caenorhabditis elegans"),
                     chicken     = c( db = "org.Gg.eg.db", gmt = "Gallus_gallus", msigdb = "Gallus gallus"),
                     chimp       = c( db = "org.Pt.eg.db", gmt = NA, msigdb = NA),
                     coelicolor  = c( db = "org.Sco.eg.db", gmt = NA, msigdb = NA), 
                     ecolik12    = c( db = "org.EcK12.eg.db", gmt = NA, msigdb = NA),
                     ecsakai     = c( db = "org.EcSakai.eg.db", gmt = NA, msigdb = NA),
                     fly         = c( db = "org.Dm.eg.db", gmt = NA, msigdb = NA),
                     gondii      = c( db = "org.Tgondii.eg.db", gmt = NA, msigdb = NA),
                     human       = c( db = "org.Hs.eg.db", gmt = "Homo_sapiens", msigdb = "Homo sapiens"),
                     malaria     = c( db = "org.Pf.plasmo.db", gmt = NA, msigdb = NA),
                     mouse       = c( db = "org.Mm.eg.db", gmt = "Mus_musculus", msigdb = "Mus musculus"),
                     pig         = c( db = "org.Ss.eg.db", gmt = "Sus_scrofa", msigdb = "Sus scrofa"),
                     rat         = c( db = "org.Rn.eg.db", gmt = "Rattus_norvegicus", msigdb = "Rattus norvegicus"),
                     rhesus      = c( db = "org.Mmu.eg.db", gmt = NA, msigdb = NA),
                     xenopus     = c( db = "org.Xl.eg.db", gmt = NA, msigdb = NA),
                     yeast       = c( db = "org.Sc.sgd.db", gmt = "Saccharomyces_cerevisiae", msigdb = "Saccharomyces cerevisiae"),
                     zebrafish   = c( db = "org.Dr.eg.db", gmt = "Danio_rerio", msigdb = "Danio_rerio")
    )
    return(annoDb)
  }
}

download_gmt <- function(species, wikipathways_folder, date = NULL){
  download <- TRUE
  species_name <- get_db(species)['gmt']
  date <- lubridate::ymd(date)
  assertthat::assert_that(assertthat::is.dir(wikipathways_folder), assertthat::is.writeable(wikipathways_folder))
  # -INF if no files correspond
  latest_species_file <- fs::dir_ls(wikipathways_folder, type = 'file') %>% 
    fs::path_file() %>%
    stringr::str_subset(species_name) %>% 
    stringr::str_extract("\\d{8}") %>% 
    lubridate::ymd() %>% 
    max()
  
  available_date <- xml2::read_html("http://data.wikipathways.org/") %>% 
    rvest::html_node("table") %>% 
    rvest::html_table() %>% 
    dplyr::filter(stringr::str_detect(Filename, '[0-9]+')) %>%
    dplyr::transmute(Date = lubridate::ymd(Filename)) %>%
    dplyr::pull() 
  
  if(!rlang::is_empty(date)) assertthat::assert_that(assertthat::is.date(date),is_in(date,available_date))
  
  latest_update_site <- max(available_date)
  
  date_char <- {if(rlang::is_empty(date)) latest_update_site else date} %>% stringr::str_remove_all("-")
  
  if (latest_species_file == latest_update_site){
    message("File on disk is already up to date")
    download <- FALSE
  } else if(latest_species_file < latest_update_site){
    file_to_download <- xml2::read_html(glue::glue("http://data.wikipathways.org/{date_char}/gmt/")) %>% 
      rvest::html_node("table") %>% 
      rvest::html_table() %>%  
      dplyr::select(Filename) %>% 
      dplyr::filter(stringr::str_detect(Filename, species_name)) %>% 
      dplyr::pull(Filename)
  }

  if(download){
    r <- httr::GET(glue::glue("http://data.wikipathways.org/{date_char}/gmt/{file_to_download}"), httr::write_disk(fs::path(wikipathways_folder,file_to_download)))
    httr::warn_for_status(r,"dowload gmt file")
    httr::message_for_status(r)
  }
} 

# "/home/aurel/Nextcloud/Documents/Cours/M2_GENIOMHE/Stage/RESSOURCES/gmt_file/wikipathways-20200310-gmt-Homo_sapiens.gmt"
load_wiki_pathway <- function(species, gmt_folder, date = NULL){
  species_name <- get_db(species)['gmt']
  
  path_species <- fs::dir_ls(gmt_folder, type = 'file') %>% 
    fs::path_file() %>% 
    stringr::str_subset(species_name) %>% 
    tibble::enframe(value = "path") %>% 
    dplyr::mutate(Date = stringr::str_extract(path,"\\d{8}")  %>% 
                    lubridate::ymd()) %>%
    dplyr::select(path, Date) %>%
    purrr::when(is.null(date) ~ dplyr::filter(.,Date == max(Date)),
                ~ dplyr::filter(.,Date == date)) %>%
    dplyr::pull(path)
  if(length(path_species) == 0){
    stop("The requested file : '", stringr::str_c("wikipathways-",as.character(date) %>% stringr::str_remove_all("-"),"-gmt-",species_name,".gmt"), "' is not found")
  } else {
    wp2gene <- clusterProfiler::read.gmt(fs::path(gmt_folder, path_species))
  }
  
  
  # BREAK changes in clusterProfiler::read.gmt : https://github.com/YuLab-SMU/clusterProfiler/blob/master/NEWS.md
  # Before 3.15.2 var name was ont, after it's term
  if (any(names(wp2gene) == "ont")){
    wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
  } else if (any(names(wp2gene) == "term")){
    wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
  }
  wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
  wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME
  return(list("TERM2GENE" = wpid2gene, 
              "TERM2NAME" = wpid2name))
}

load_msigdb <- function(species, category){
  species_name <- get_db(species)['msigdb']
  ### Signatures Moléculaires: MSigDB
  #Avec un fichier de ref d'une signature: 
  #H: hallmark gene sets
  #C1: positional gene sets
  #C2: curated gene sets
  #C3: motif gene sets
  #C4: computational gene sets
  #C5: GO gene sets
  #C6: oncogenic signatures
  #C7: immunologic signatures
  #msigdbr_show_species() #montre toutes les espèces possibles
  m_t2g <- msigdbr::msigdbr(species = species_name, category = category) %>% dplyr::select(gs_name, entrez_gene)
  return(m_t2g)
}

load_cell_markers <- function(species){
  cell_markers <- suppressMessages(
    glue::glue('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/{stringr::str_to_title(species)}_cell_markers.txt') %>%
      vroom::vroom() %>%
      tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
      dplyr::select(cellMarker, geneID) %>%
      dplyr::mutate(geneID = strsplit(geneID, ', ')))
  return(cell_markers)
} 

test_if_done <- function(test, DE_type, output_path_def){
  path_comp <- output_path_def %>% fs::path_dir() %>% fs::path("Comparisons")
  comp <- fs::dir_ls(path_comp) %>% 
    fs::path_filter(regexp = test) %>%
    tibble::as_tibble_col(column_name = "path") %>%
    purrr::when(nrow(.) == 0 ~ tibble::tibble(method = "comparisons", test = test, data = FALSE, Plot = FALSE),
                TRUE ~ dplyr::mutate(.,dir = fs::is_dir(path), file = fs::is_file(path), name = fs::path_file(path), ext = fs::path_ext(path)) %>%
                  dplyr::group_by(ext) %>%
                  dplyr::count() %$%
                  tibble::tibble(method = "comparisons", test = test, data = ifelse(n == 1 & ext == "rds", TRUE, FALSE), Plot = ifelse(n >= 1 & ext == "png", TRUE, FALSE)) %>%
                  dplyr::group_by(method,test) %>%
                  dplyr::summarise(data = any(data), Plot = any(Plot), .groups = "drop"))
    
  
  out <- tibble::tibble(method = c('enrichment', 'GSEA'), test = test)
  
  list_file <- fs::dir_ls(output_path_def) %>% 
    tibble::as_tibble_col(column_name = "path") %>%
    dplyr::mutate(dir = fs::is_dir(path), file = fs::is_file(path), name = fs::path_file(path), ext = fs::path_ext(path)) 
  
  suppressMessages(out <- list_file %>% 
                     dplyr::filter(ext == "rds") %>%
                     dplyr::select(name) %>%
                     purrr::when(nrow(.) == 0 ~ tibble::tibble(method = c("enrichment", 'GSEA'), data = FALSE),
                                 TRUE ~ dplyr::transmute(.,name = fs::path_ext_remove(name) %>% stringr::str_extract('enrichment|gsea')) %>%
                                   tidyr::unnest_wider(col = "name") %>% 
                                   dplyr::rename(method = ...1) %>%
                                   dplyr::mutate(data = TRUE) %>%
                                   dplyr::bind_rows(., dplyr::mutate(.,method = stringr::str_to_upper(method))) ) %>%
                     dplyr::right_join(out, by = c('method')))
  n_files <- out %>% 
    purrr::when(!all(.$data) ~ c("enrichment" = NA , 'GSEA' = NA),
                TRUE ~ list_file %>% 
                  dplyr::filter(stringr::str_detect(name,"plot")) %>% 
                  dplyr::pull(path) %>% 
                  fs::path("n_files.txt") %>% 
                  as.vector() %>% 
                  magrittr::set_names(stringr::str_extract(.,"Enrichment|GSEA") %>% stringr::str_replace("Enrich", "enrich")) %>%
                  purrr::map_dbl(~purrr::possibly(readr::read_lines, otherwise = NA)(.) %>% as.double()) )
  
  suppressMessages(out <- list_file %>%
                     dplyr::filter(dir) %>%
                     dplyr::pull(path) %>%
                     fs::dir_ls(recurse = TRUE) %>%
                     tibble::as_tibble_col("path") %>%
                     dplyr::mutate(ext = fs::path_ext(path), filename = fs::path_file(path) %>% fs::path_ext_remove(), dir_1 = fs::path_dir(path) %>% fs::path_file(), dir_2 = fs::path_dir(path) %>% fs::path_dir() %>% fs::path_file()) %>%
                     dplyr::group_by(dir_1, dir_2, .drop = FALSE) %>% 
                     dplyr::count() %>% dplyr::filter((dir_2 == "Enrichment_plot" | dir_2 == "GSEA_plot") & dir_1 != "Pathways") %>%
                     purrr::when(nrow(.) == 0 ~ tibble::tibble(method = c("enrichment", 'GSEA'), Plot = FALSE),
                                 TRUE ~ dplyr::mutate(.,n_files = dplyr::if_else(dir_2 == "Enrichment_plot", n_files['enrichment'], n_files['GSEA'])) %>%
                                   dplyr::mutate(Plot = dplyr::if_else(n == n_files, TRUE, FALSE), method = stringr::str_split(dir_2,'_')) %>% 
                                   tidyr::unnest_wider(method) %>%
                                   dplyr::rename(method = ...1) %>%
                                   dplyr::select(method,Plot) %>%
                                   dplyr::group_by(method) %>%
                                   dplyr::summarise(Plot = any(Plot), .groups = "drop") %>%
                                   dplyr::bind_rows(., dplyr::mutate(.,method = stringr::str_to_lower(method)))) %>%
                     dplyr::right_join(out, by = c('method')) %>%
                     tidyr::replace_na(list(data = FALSE, test = "", method = "", Plot = FALSE)))
  
  return(dplyr::bind_rows(out,comp))
}

result_message_summary <- function(data, path){
  purrr::map(data, ~ purrr::pluck(.,"result")) -> out
  purrr::map(data, ~ purrr::pluck(.,"messages")) %>% purrr::compact() -> mess
  purrr::map(data, ~ purrr::pluck(.,"warnings")) %>% purrr::compact() -> warn
  
  unlist(mess) %>%
    tibble::enframe(name = 'method', value = "mess") %>%
    dplyr::mutate(no_term_enrich = ifelse(stringr::str_detect(mess,"no term enriched under specific pvalueCutoff"),TRUE,FALSE),
                  no_set_size = ifelse(stringr::str_detect(mess,"No gene set have size "), TRUE, FALSE),
                  no_map = ifelse(stringr::str_detect(mess,"No gene can be mapped...."), TRUE, FALSE)) %>%
    dplyr::group_by_at(dplyr::vars(-mess)) %>% 
    dplyr::count() %>% 
    dplyr::mutate(message = dplyr::case_when(no_term_enrich ~ paste0("\t\t\t",n, " ", method, " anlysis gave no term enriched under specific pvalueCutoff"),
                                             no_set_size ~ paste0("\t\t\tGene set do not have the required size ",n, " times in ", method),
                                             no_map ~ paste0("\t\t\tNo genes can be mapped ",n, " times in ", method, " returning NULL"))) %>% 
    dplyr::pull(message) %>% 
    paste0(collapse = '\n') %>% 
    stringr::str_remove_all("NA|\nNA") -> to_print

  if(to_print != "") message(to_print)
  
  # use case_when for different message
  unlist(warn) %>% 
    tibble::enframe(name = 'method', value = "warn") %>%
    dplyr::mutate(pval_not_calc = ifelse(stringr::str_detect(warn,"pathways for which P-values were not calculated properly due to unbalanced gene-level statistic values"),TRUE,FALSE),
                  n_perm_used = ifelse(stringr::str_detect(warn,"We do not recommend using nPerm parameter"),TRUE,FALSE),
                  ties =  ifelse(stringr::str_detect(warn,"The order of those tied genes will be arbitrary, which may produce unexpected results."),TRUE,FALSE),
                  duplicates = ifelse(stringr::str_detect(warn,"There are duplicate gene names, fgsea may produce unexpected results."),TRUE,FALSE)) %>%
    dplyr::group_by_at(dplyr::vars(-warn)) %>% 
    dplyr::count() %>% 
    dplyr::mutate(message = dplyr::case_when(pval_not_calc ~ paste0("\t\t\t",n, " P-values were not calculated properly in ", method),
                                             n_perm_used ~ paste0("\t\t\t","nPerm parameter is used ",n, " times in ", method, " which is not recommended"),
                                             ties ~ paste0("\t\t\t","There were ties ",n, " times in ", method, " which may produce unexpected results"),
                                             duplicates ~ paste0("\t\t\t","There were duplicates gene names ",n, " times in ", method, " which may produce unexpected results"))) %>% 
    dplyr::pull(message) %>% 
    paste0(collapse = '\n') %>% 
    stringr::str_remove_all("NA|\nNA") -> to_print 
  
  if(to_print != "") message(to_print)
  
  n_files <- purrr::map(out, names) %>% unlist %>% stringr::str_subset("Ontology|Pathway|Signatures|Cell") %>% length()
  readr::write_lines(n_files , path = fs::path(path,"n_files.txt"))
  
  return(out)
}

viewPathway_igraph <- function(pathName, organism = "human", readable = TRUE, foldChange = NULL, keyType = "ENTREZID"){
  
  
  # convertion to the names that graphite::pathways understands
  org2org <- list(arabidopsis="athaliana",
                  bovine="btaurus",
                  canine="cfamiliaris",
                  chicken="ggallus",
                  ecolik12="ecoli",
                  fly="dmelanogaster",
                  human="hsapiens",
                  mouse="mmusculus",
                  pig="sscrofa",
                  rat="rnorvegicus",
                  celegans="celegans",
                  xenopus="xlaevis",
                  yeast="scerevisiae",
                  zebrafish="drerio")
  
  if(!(organism %in% names(org2org))){
    cat(paste(c("the list of supported organisms:",names(org2org)), collapse='\n'))
    stop(sprintf("organism %s is not supported", organism))
  }
  
  p <- graphite::pathways(org2org[[organism]], 'reactome')[[pathName]]
  
  if (readable) {
    p <- purrr::quietly(graphite::convertIdentifiers)(p, "symbol")[['result']]
    if (!is.null(foldChange)) {
      stopifnot(!any(duplicated(names(foldChange)))) # can't have two value for one gene
      OrgDb <- get_db(organism)['db']
      names(foldChange) <- DOSE::EXTID2NAME(OrgDb, names(foldChange), keyType)
      
    }
  } else {
    if (!is.null(foldChange)) {
      p <- graphite::convertIdentifiers(p, "entrez")
    }
  }
  
  gg <- graphite::pathwayGraph(p) %>%
    igraph::igraph.from.graphNEL() %>%
    igraph::as.undirected() %>%
    igraph::set_vertex_attr("size", value = 8) %>%
    igraph::set_vertex_attr("color", value = "#B3B3B3") %>%
    igraph::set_vertex_attr("label", value = igraph::vertex_attr(.,"name")) %>%
    igraph::set_edge_attr("width", value = 1) %>%
    # igraph::set_edge_attr("color", value = "#8DA0CB") %>%
    igraph::set_vertex_attr("name", value = sub("[^:]+:", "",  igraph::vertex_attr(.,"name")))
  
  if (!is.null(foldChange)) {
    
    fc_color <- foldChange[igraph::V(gg)$name] %>%
      scales::cscale(scales::seq_gradient_pal("red", "blue"), na.value = "#E5C494")
    gg <- igraph::set_vertex_attr(gg,"color", value = fc_color) %>%
      igraph::set_vertex_attr("color_value", value = foldChange[igraph::V(gg)$name])
    
  }
  
  return(gg)
}

plot_igraph_interactive <- function(gg, title, seed = 1337, params){
  net <- visNetwork::toVisNetworkData(gg, idToLabel = TRUE) %$%
    visNetwork::visNetwork(nodes = nodes, edges = edges, main = title, submain = paste0("Cluster ", params$tested_cluster, " vs ", params$control_cluster, " for min.pct ", params$min.pct), width = 1800, height = 800) %>%
    visNetwork::visIgraphLayout(layout = 'layout_with_lgl',type = 'full', randomSeed = seed) %>%
    visNetwork::visEdges(color = list(color = "rgba(200,200,200,0)", highlight = "rgba(140,140,140,0.4)", hover = "blue")) %>%
    visNetwork::visOptions(selectedBy = list(variable = "label", hideColor = "rgba(200,200,200,0)", highlight = TRUE, main = "Select by gene", 
                                             style = 'width: 200px; height: 26px;
                                                      background: #f8f8f8;
                                                      color: dark;
                                                      border:none;
                                                      outline:none;'), 
                           highlightNearest = list(enabled = TRUE, hideColor = "rgba(200,200,200,0)",
                                                   degree = 1, algorithm = "all", labelOnly = FALSE))
  return(net)
}

plot_igraph_static <- function(gg, title, params, layout = 'kk'){
  p <- ggraph::ggraph(gg, layout=layout) +
    ggraph::geom_edge_link(alpha=.8, colour='darkgrey') +
    ggraph::geom_node_point(ggplot2::aes_(color=~as.numeric(as.character(color_value)), size=~size)) +
    ggplot2::scale_color_continuous(low="red", high="blue", name = "fold change", na.value = "#E5C494") +
    ggraph::geom_node_text(ggplot2::aes_(label=~name), repel=TRUE) +
    ggplot2::scale_size(guide = FALSE) + 
    ggplot2::theme_void() + 
    ggplot2::labs(title = title, subtitle = paste0("Cluster ", params$tested_cluster, " vs ", params$control_cluster, " for min.pct ", params$min.pct))
  return(p)
}

assing_out_values <- function(data,out,slot){
  data %>%
    purrr::map( ~ purrr::pluck(.,slot)) %>%
    purrr::compact() %>%
    purrr::iwalk(function(value,name) out[[slot]][[name]] <<- value)
  return(out)
}
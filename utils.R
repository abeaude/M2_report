# Set options
pbapply::pboptions(type = "none")

# Re exports
is_in <- magrittr::is_in
`%<>%` <- magrittr::`%<>%`
`%>%` <- magrittr::`%>%`
`%$%` <- magrittr::`%$%`
`%T>%` <- magrittr::`%T>%`
`%||%` <- rlang::`%||%`
is_named <- rlang::is_named

# Custom assertion
is_seurat <- function(x) {
  is(x,"Seurat")
}

are_seurat <- function(x) {
  all(purrr::map_lgl(x,is_seurat))
}

are_in <- function(x,table){
  all(is_in(x,table))
}

has_metadata <- function(sobj, metadata){
  assertthat::assert_that(is_seurat(sobj))
  metadata %in% colnames(sobj@meta.data)
}

has_cell <- function(sobj,cell){
  assertthat::assert_that(is_seurat(sobj))
  cell %in% Seurat::Cells(sobj)
}

has_dimred <- function(sobj,dimred){
  assertthat::assert_that(is_seurat(sobj))
  all(purrr::map_lgl(dimred,~any(stringr::str_detect(.,stringr::regex(Seurat::Reductions(sobj), ignore_case = TRUE)))))
}

has_assay <- function(sobj,assay){
  assay %in% names(sobj@assays)
}

assertthat::on_failure(has_assay) <- function(call, env) {
  paste0(deparse(call$sobj), " does not have an assay named ", deparse(call$assay))
}

assertthat::on_failure(is_seurat) <- function(call, env) {
  paste0(deparse(call$x), " is not a Seurat object")
}

assertthat::on_failure(are_seurat) <- function(call, env) {
  list_sobj <- eval(call$x,env)
  is_seurat_obj <- purrr::map_lgl(list_sobj,is_seurat)
  wrong_elements <- which(isFALSE(is_seurat_obj)) %>% paste(collapse = ",")
  paste0("Elements ",wrong_elements, " of ", deparse(call$x), " are not a Seurat object")
}

assertthat::on_failure(has_metadata) <- function(call, env) {
  paste0(deparse(call$metadata), " is not in the metadata of the Seurat object ", deparse(call$sobj))
}

assertthat::on_failure(has_cell) <- function(call, env) {
  paste0(deparse(call$cell), " is not a cell of the Seurat object ", deparse(call$sobj))
}

assertthat::on_failure(has_dimred) <- function(call, env) {
  dimred <- eval(call$dimred,env)
  sobj_dimred <- eval(call$sobj,env)
  wrong_dimred <- paste0(purrr::map(purrr::set_names(dimred,dimred),~stringr::str_detect(.,stringr::regex(Seurat::Reductions(sobj), ignore_case = TRUE))) %>% purrr::modify(any) %>% purrr::keep(isFALSE) %>% names(),collapse = ",")
  paste0("Seurat object `", deparse(call$sobj),"` does not have reductions dimensions : ", wrong_dimred)
}

assertthat::on_failure(is_named) <- function(call, env) {
  paste0(deparse(call$x), " is not a named ", class(eval(call$x,env)))
}

assertthat::on_failure(is_in) <- function(call, env) {
  paste0(deparse(call$x), " is not in ", deparse(call$table))
}

assertthat::on_failure(are_in) <- function(call, env) {
  x <- eval(call$x,env)
  table <- eval(call$table,env)
  are_in_res <- purrr::map2_lgl(x,table,is_in) %>% purrr::set_names(x)
  wrong_elements <- are_in_res[are_in_res == FALSE] %>% names() %>% paste(collapse = ",")
  paste0("Elements ", wrong_elements, " of ", deparse(call$x), " are not in ", deparse(call$table) %>% stringr::str_c(collapse = "") %>% stringr::str_squish())
}

# Source a folder
sourceFolder <- function(dir){
  fs::dir_walk(path = dir, fun = source, type = "file") 
}

# Conversion
seurat_to_anndata <- function(sobj){
  reticulate::use_condaenv("scanpy")
  sc <- reticulate::import("scanpy")
  
  adata_seurat <- sc$AnnData(
    X   = Matrix::t(Seurat::GetAssayData(sobj)),
    obs = sobj[[]],
    var = Seurat::GetAssay(sobj)[[]]
  )
  
  dimred <- Seurat::Reductions(sobj)
  
  adata_seurat$obsm <- purrr::map(dimred, ~ Seurat::Embeddings(sobj,.x)) %>% purrr::set_names(dimred) %>%
    reticulate::dict()
  adata_seurat
}

seurat_is_normalized <- function(sobj){
  sobj@commands %>% names() %>% stringr::str_detect("SCTransform|NormalizeData") %>% any()
}

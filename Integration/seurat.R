seurat_SCT_integration <- function(){
  
  assay <- rep(assay,length(sobj_list))
  # Use VariableFeatures if previously computed in the object
  # sobj_integrated@assays$RNA@var.features
  
  message("Integration preparation")
  sobj.features <- Seurat::SelectIntegrationFeatures(object.list = sobj_list, nfeatures = 3000, assay = assay, verbose = TRUE)
  sobj_list <- Seurat::PrepSCTIntegration(object.list = sobj_list, assay = assay, anchor.features = sobj.features, verbose = TRUE)
  sobj.anchors <- Seurat::FindIntegrationAnchors(object.list = sobj_list, assay = assay_vec,reference = reference, normalization.method = "SCT", anchor.features = sobj.features, verbose = TRUE, scale = FALSE, reduction = 'cca',l2.norm = TRUE)
  
  sobj_integrated <- Seurat::IntegrateData(anchorset = sobj.anchors, normalization.method = "SCT", verbose = TRUE)
  Seurat::DefaultAssay(sobj_integrated) <- "integrated"
  return(sobj_integrated)
}

seurat_standard_integration <- function(sobj_list, n_features, assay, verbose){
  sobj_list <- purrr::map(sobj_list, ~Seurat::FindVariableFeatures(object = .x, selection.method = "vst", nfeatures = n_features, verbose = verbose))
  
  sobj.anchors <- Seurat::FindIntegrationAnchors(object.list = sobj_list, assay = assay, reference = NULL, normalization.method = "LogNormalize", anchor.features = n_features, verbose = TRUE, scale = TRUE, reduction = 'cca', l2.norm = TRUE, dims = dims, k.anchor = 5,k.filter = 200,k.score = 30,max.features = 200, nn.method = "rann",)
  
  sobj_integrated <- Seurat::IntegrateData(anchorset = sobj.anchors, normalization.method = "LogNormalize", dims = dims, verbose = TRUE)
  
  Seurat::DefaultAssay(sobj_integrated) <- "integrated"
  return(sobj_integrated)
}
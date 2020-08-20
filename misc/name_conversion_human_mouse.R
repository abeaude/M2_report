convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  # humanx <- unique(genesV2[, 2])
  
  genesV2
  # Print the first 6 genes found to the screen
  # return(humanx)
}

convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

genes <- convertMouseGeneList(mouse.genes)
# saveRDS(genes, file = "/mnt/Data/RESSOURCES/GENELISTS/homo_sapiens_stress_core_genes.rds")

TCELL.markers <- c("CD3G", "CD3D", "CD3E") ## T-cells c("Cd3g", "Cd3d", "Lat", "Cd3e", "Skap1", "Il7r", "Foxp3")
names(TCELL.markers) <- rep("T Cell", length(TCELL.markers))
TH.markers <- c("CD4")
names(TH.markers) <- rep("T Helpers and regulatory",length(TH.markers))
TEFF.markers <- c("CD8A","GZMB","GZMA","FASLG")
names(TEFF.markers) <- rep("Effector T cells",length(TEFF.markers))
TREG.markers <- c("FOXP3","IL2RA","IL7R")
names(TREG.markers) <- rep("Treg",length(TREG.markers))
SCHWANN.markers <- c("SOX10", "PLP1", "CADM4", "NGFR", "POU3F1","ANK3") ## Schwann
names(SCHWANN.markers) <- rep("Schwann Cell", length(SCHWANN.markers))
FIBRO.markers <- c("PDGFRA", "PDGFRB", "DCN", "MGP", "SERPING1","CXCL14") ## Fibroblasts
names(FIBRO.markers) <- rep("Fibroblast", length(FIBRO.markers))
ENDO.markers <- c("PECAM1", "CD34","EGFL7", "EMCN") ## Endothelial
names(ENDO.markers) <- rep("Endothelial cell", length(ENDO.markers))
MUSCLE.markers <- c("RGS5", "SERPINA1", "HIGD1B") ## Smooth muscle
names(MUSCLE.markers) <- rep("Smooth muscle cell", length(MUSCLE.markers))
IMMUNOCOM.markers <- c("PTPRC") ## Common to all immune
names(IMMUNOCOM.markers) <- rep("Immunity (common)", length(IMMUNOCOM.markers))
MACRO.markers <- c( "ADGRE1", "CD68", "LYZ") ## Macrophages c("Adgre1", "P2ry6", "C3ar1", "Maf", "Ms4a7", "Cd68",Lyz2)
names(MACRO.markers) <- rep("Macrophage", length(MACRO.markers))
MACRO_CMH.markers <- c( "C1QA", "C1QB", "C1QC", "TREM2", "APOE")
names(MACRO_CMH.markers) <- rep("Macrophages Complement +",rep(length(MACRO_CMH.markers)))
MACRO_plac8.markers <- c( "PLAC8")
names(MACRO_plac8.markers) <- rep("Macrophages Plac8+",length(MACRO_plac8.markers))
MACRO_cd209a.markers <- c( "CD209A")
names(MACRO_cd209a.markers) <- rep("Macrophages Cd209a+",length(MACRO_plac8.markers))
NK.markers <- c("XCL1", "KLRK1", "NCR1", "KLRB1","NKG7") ## NK cells
names(NK.markers) <- rep("Natural Killer", length(NK.markers))
BCELL.markers <- c("CD79A", "LY6D", "CD79B") ## B cells
names(BCELL.markers) <- rep("B cell", length(BCELL.markers))
DC.markers <- c("XCR1") ## DC cells
names(DC.markers) <- rep("Dendritic cell", length(DC.markers))
pDC.markers <- c("CCR9","CD33") ## DC cells
names(pDC.markers) <- rep("Dendritic cell", length(pDC.markers))
NEUTROPHIL.markers <- c( "S100A9", "S100A8") #neutrophils
names(NEUTROPHIL.markers) <- rep("Neutrophils", length(NEUTROPHIL.markers))

# scriptTemplate.R
#
# Purpose: To parse data from the Human Protein Atlas TSV and XML files
# Version: 1.0
# Date: 2019-02-02
# Author: Nada Elnour
# License:
#
# Input:
# Output:
# Dependencies: XML 3.98-1.16; stringr 1.3.1; tibble 2.0.1; biomaRt 2.38.0; grImport 0.9-1.1; rsvg 1.3
# ==============================================================================

# NO SIDE EFFECTS:
# This script can be safely source()'d to define the functions it contains and
# install.packages()/run library() as required.
# All other code will not be executed unless this is done interactively.


# ====  PARAMETERS  ============================================================
# Define and explain all parameters. No "magic numbers" in your code below.


#
#      This script uses guard-blocks that prevent execution of
#      code that should not be executed when the entire script
#      is sourced. Thus it can be source()'d to load its functions,
#      or executed interactively.
#
if (FALSE) { # <--- Keep this guard-block only if your script needs to be
  #      source()'d without side effects, since loading parameters
  #      changes the global workspace. If side-effects are oK
  #      remove this guard block, and the others in this script as
  #      required.
  
}


# ====  PACKAGES  ==============================================================
# Load all required packages.

if (requireNamespace("tibble", quietly=TRUE)) {
  library("tibble")
} else {
  install.packages("tibble")
  library(tibble)
}

if (requireNamespace("stringr", quietly=TRUE)) {
  library("stringr")
} else {
  install.packages("stringr")
  library(stringr)
}

if (! requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("biomaRt", version = "3.8")
}
library(biomaRt)

# ====  FUNCTIONS  =============================================================

mapColors <- function(cell_picture){
  # Purpose: Maps color HEX codes in <cell_picture> to subcelular loci
  # Parameters:
  #     cell_picture: A Picture class object of the cell
  # Value:
  #     result: A dataframe mapping the specifying which HEX color code
  #             stands for which subcellular component 
  
  allPaths <- explodePaths(cell_picture)
  colorCodes <- data.frame(RGB = character())
  
  for (paths in allPaths@paths){
    colorCodes <- rbind(colorCodes, data.frame(RGB=paths@rgb))
  }
  
  colorCodes <- unique(colorCodes)
  
  colorCodes <- cbind(colorCodes, data.frame(Loci = c("Cytosol",
                                                      "Cytosol",
                                                      "Actin",
                                                      "Intermediate filaments",
                                                      "Outline",
                                                      "Actin",
                                                      "Actin",
                                                      "Actin",
                                                      "Actin",
                                                      # "Nuclear membrane",
                                                      "Microtubule organizing center",
                                                      "Centrosome",
                                                      "Microtubules",
                                                      "Microtubules",
                                                      "Lipid droplets",
                                                      "Lysosomes",
                                                      "Peroxisomes",
                                                      "Endosomes",
                                                      "Endoplasmic reticulum",
                                                      "Golgi apparatus",
                                                      "Nucleoplasm",
                                                      "Nucleus",
                                                      "Nucleoli",
                                                      "Nuclear speckles",
                                                      "Nuclear bodies",
                                                      "Nucleoli fibrillar center",
                                                      "Rods and rings",
                                                      "Mitochondria",
                                                      "Mitochondria",
                                                      "Plasma membrane",
                                                      "Plasma membrane"
  )))
  
  save(colorCodes, file = "./inst/extdata/colorCodes.RData")
}

parseHPAData <- function(filepath, reliability = "enhanced") {
  # Purpose: Parse Human Protein Atlas TSV file into a dataframe of mapped HGNC symbols 
  # Parameters:
  #     filepath: full path reference to TSV file
  #     reliability: reliability score category for filtration; one of: "enhanced", "approved", "uncertain", or "supported".
  # Value:
  #     result: A dataframe wth gene names, ENSEMBL IDs, NCBI Refseq IDs, localization information, GO IDs, HGNC symbols, and chromosome loci
  
  data <- read.delim(filepath)
  
  if (reliability == "enhanced"){
    data <- data[data$Reliability == "Enhanced", ]
    toDelete <- c("Reliability", "Supported", "Approved", "Uncertain")
    data <- data[, !(colnames(data) %in% toDelete), drop=FALSE]
  }
  else if (reliability == "approved"){
    data <- data[data$Reliability == "Approved", ]
    toDelete <- c("Reliability", "Supported", "Enhanced", "Uncertain")
    data <- data[, !(colnames(data) %in% toDelete), drop=FALSE]
  }
  else if (reliability == "supported"){
    data <- data[data$Reliability == "Supported", ]
    toDelete <- c("Reliability", "Enhanced", "Approved", "Uncertain")
    data <- data[, !(colnames(data) %in% toDelete), drop=FALSE]
  }
  else if (reliability == "uncertain"){
    data <- data[data$Reliability == "Uncertain", ]
    toDelete <- c("Reliability", "Supported", "Approved", "Enhanced")
    data <- data[, !(colnames(data) %in% toDelete), drop=FALSE]
  }
  toAdd <- getHGNCRefs(data$Gene)
  
  # check if there are symbols in toAdd not in data and remove them
  if (length(data$Gene) > length(toAdd$ensembl_gene_id)){
    diff <- data$Gene[!(data$Gene %in% toAdd$ensembl_gene_id)]
    data <- data[!(data$Gene %in%  diff), ]
  } else if (length(data$Gene) < length(toAdd$ensembl_gene_id)) {
    diff <- toAdd$Approved.symbol[!(toAdd$Approved.symbol %in% data$Gene.name | 
                                      toAdd$Previous.symbols %in% data$Gene.name | 
                                      toAdd$Synonyms %in% data$Gene.name)]
    toAdd <- toAdd[!(toAdd$Approved.symbol %in%  diff), ]
  }
  
  data <- merge(data, toAdd, by.x="Gene", by.y = "ensembl_gene_id")
  data <- data[!duplicated(data$hgnc_symbol), ]
  rownames(data) <- data$hgnc_symbol
  return(data)
}

getHGNCRefs <- function(gene_id){
  # Purpose: Searches hsapien dataset in ENSEMBL biomaRt for gene names, ENSEMBL IDs, HGNC symbols, chromsome loci, 
  #          and peptide Refseq IDs of the genes referred to by gene_id
  # Parameters:
  #     gene_id: list of ENSEMBL ID strings
  # Value:
  #     result: A dataframe wth gene names, ENSEMBL IDs, HGNC symbols, chromsome loci, and peptide Refseq IDs for the genes in gene_id
  ensembl <- useMart(biomart = "ensembl")
  human <- searchDatasets(mart = ensembl, pattern = "hsapiens")
  
  myMart <- useMart("ensembl", human$dataset)
  att <- c("ensembl_gene_id", "hgnc_symbol", "external_gene_name", "chromosome_name", "refseq_peptide")
  dataRefs <- getBM(attributes = att, mart = myMart,
                    filters = "ensembl_gene_id",
                    values = gene_id)
  
  dataRefs <- dataRefs[!duplicated(dataRefs$ensembl_gene_id), ]
  return(dataRefs)
}

whereIs <- function(hgnc_symbol, HPASet){
  # Purpose: Search and colour subcellular locations where protein of interest is annotated to localize
  # Parameters:
  #     hgnc_symbol: HGNC symbol of gene of protein of interest
  #     HPASet: Human Protein Atlas dataset
  # Value:
  #     result: an image of the cel highlighted for subcellular components where protein of gene interest is annotated to localize
  if (requireNamespace("XML", quietly=TRUE)) {
    library("XML")
  } else {
    install.packages("XML")
    library(XML)
  }
  
  if (requireNamespace("grid", quietly=TRUE)) {
    library("grid")
  } else {
    install.packages("grid")
    library(grid)
  }
  
  if (requireNamespace("grImport", quietly=TRUE)) {
    library("grImport")
  } else {
    install.packages("grImport")
    library(grImport)
  }
  
  if (requireNamespace("rsvg", quietly=TRUE)) {
    library("rsvg")
  } else {
    install.packages("rsvg")
    library(rsvg)
  }
  
  # Reference cell: provided image of the cell was obtained from Human Protein Atlas download page
  # https://www.proteinatlas.org/images_static/cell.svg
  rsvg_ps("./inst/extdata/cell.svg", "./inst/extdata/cell.ps")
  PostScriptTrace("./inst/extdata/cell.ps", outfilename = "./inst/extdata/cell.xml")
  cell <- readPicture("./inst/extdata/cell.xml")
  
  geneLoci <- unlist(strsplit(as.character(test[hgnc_symbol,]$Enhanced),";"))
  
  load("./inst/extdata/colorCodes.RData")
  
  toChange <- colorCodes$RGB[colorCodes$Loci %in% geneLoci]
  
  for (idx in 1:length(cell@paths)){
    if (cell@paths[idx]$path@rgb %in% toChange){
      cell@paths[idx]$path@rgb <-  stringr::str_replace(cell[loci]@paths$path@rgb, ".{7}", "#010202")
    }
  }
  
  #plot the good stuff
  grid.picture(cell)
}

# ====  PROCESS  ===============================================================
# Enter the step-by-step process of your project here. Strive to write your
# code so that you can simply run this entire block and re-create all
# intermediate results.
if (FALSE) {
  
  filepath = "../data/subcellular_location.tsv"
  enhanced <- parseHPAData(filepath)
  
  # What are the most commonly annotated subcellular localization sites?
  par(las=2)
  par(mar=c(3,15,0,1))
  barplot(table(test$Enhanced), horiz = TRUE, cex.names = 0.5, cex.axis = 0.8)
  
  # Visualize a gene's protein localization information
  whereIs("DVL2", enhanced)
  # which localizes to 
  unlist(strsplit(as.character(enhanced["DVL2",]$Enhanced),";"))
  
}

# ====  TESTS  =================================================================
if (FALSE) {
  # Enter your function tests here...
}


# [END]
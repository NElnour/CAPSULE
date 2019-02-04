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
# Dependencies:
#
# ToDo:
# Notes:
#
# ==============================================================================

# Keep only one of the two notices below, and then remove this line.

# WARNING: SIDE EFFECTS
# Executing this script will execute code it contains.

# NO SIDE EFFECTS:
# This script can be safely source()'d to define the functions it contains and
# install.packages()/run library() as required.
# All other code will not be executed unless this is done interactively.




# Use setwd() with discretion - it should normally not be necessary to change
# the working directory, and it is poor practice since it changes the global
# state. If you must use setwd(), save the current working directory and restore
# it when your script is done, as per the example code below.
#
# oldWD <- getwd()
# setwd("<your/project/directory>")
#
# setwd(oldWD)   # <--- move this to the end of your script

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

if (requireNamespace("seqinr", quietly=TRUE)) {
  library("seqinr")
} else {
  install.packages("seqinr")
  library(seqinr)
}
# Package information:
#  library(help = seqinr)       # basic information
#  browseVignettes("seqinr")    # available vignettes
#  data(package = "seqinr")     # available datasets

if (requireNamespace("xml2", quietly=TRUE)) {
  library("xml2")
} else {
  install.packages("xml2")
  library(xml2)
}

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
}
if (! requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}

# ====  FUNCTIONS  =============================================================

# Define functions or source external files
if (FALSE) { # <---- If your script needs to be side-effect proof, source
  #       only scripts that are themselves side-effect proof or
  #       keep this guard block in place.
  
  source("<myUtilityFunctionsScript.R>")
  
}

parseHPAXML <- function(filepath) {
  # Purpose: Parse Human Protein Atlas XML file into a dataframe of reliability scores
  #     Describe ...
  # Parameters:
  #     filepath: full path reference to XML file
  # Value:
  #     result: A dataframe wth gene names, and their reliability scores
  
  # code ...
  
  data <- read_xml("./inst/extdata/proteinatlas.xml")
  
  return(result)
}

parseHPADataL <- function(filepath, reliability = "all") {
  # Purpose: Parse Human Protein Atlas XML file into a dataframe of reliability scores
  #     Describe ...
  # Parameters:
  #     filepath: full path reference to XML file
  # Value:
  #     result: A dataframe wth gene names, and their reliability scores
  
  # code ...
  
  data <- read.delim("./inst/extdata/subcellular_location.tsv")
  data <- data[!duplicated(data$Gene.name), ]
  rownames(data) <- data$Gene.name
  
  HGNC <- getHGNCRefs()$HGNCData
  withdrawnSymbols <- getHGNCRefs()$withdrawnSymbols
  
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
  
  # select all HGNC columns whose rows match with the gene names already in data
  approvedSelection <- HGNC[HGNC$Approved.symbol %in% data$Gene.name, ]
  # if there are no matches between data gene symbols and HGNC approved symbols, check if HPA is using an HGNC synonym
  synonymousSelection <- HGNC[HGNC$Synonyms %in% data$Gene.name,]
  # check if there are sybols in data not in HGNC
  # data$Gene.name[!(data$Gene.name %in% toAdd$Approved.symbol)]
  # if there are no matches between approved and synonymous symbols, check if HPA is using a withdrawn symbol
  withdrawnSelection <- HGNC[withdrawnSymbols$idx[which(withdrawnSymbols$str_extract.dataRefs.Approved.symbol.idx.........withdrawn... %in% data$Gene.name)], ]
  # if there are no matches and the symbol is not withdrawn, check previous symbols
  previousSelection <- HGNC[HGNC$Previous.symbols %in% data$Gene.name, ]
  
  # Check if there are overlaps in the selections
  previousSelection <- previousSelection[!previousSelection$Approved.symbol %in% synonymousSelection$Approved.symbol, ]
  previousSelection <- previousSelection[!previousSelection$Approved.symbol %in% approvedSelection$Approved.symbol, ]
  previousSelection <- previousSelection[!previousSelection$Approved.symbol %in% withdrawnSelection$Approved.symbol, ]
  
  synonymousSelection <- synonymousSelection[!synonymousSelection$Approved.symbol %in% approvedSelection$Approved.symbol, ]
  synonymousSelection <- synonymousSelection[!synonymousSelection$Approved.symbol %in% withdrawnSelection$Approved.symbol, ]
  
  withdrawnSelection <- withdrawnSelection[!(withdrawnSelection$Approved.symbol %in% approvedSelection$Approved.symbol), ]
  
  toAdd <- rbind(previousSelection, synonymousSelection, withdrawnSelection, approvedSelection)
  
  # check if there are symbbbols in toAdd not in data and remove them
  diff <- toAdd$Approved.symbol[!(toAdd$Approved.symbol %in% data$Gene.name | toAdd$Previous.symbols %in% data$Gene.name | toAdd$Synonyms %in% data$Gene.name)]
  toAdd <- toAdd[!(toAdd$Approved.symbol %in%  diff), ]
  
  add_column(data, toAdd, .after = data$Gene.name)
  return(result)
}

getHGNCRefs <- function() {
  # Purpose: Parse Human Protein Atlas XML file into a dataframe of reliability scores
  #     Describe ...
  # Parameters:
  #     filepath: full path reference to XML file
  # Value:
  #     result: A dataframe wth gene names, and their reliability scores
  
  # code ...
  
  
  dataRefs <- read.delim("./inst/extdata/HGNC.txt")
  dataRefs <- dataRefs[!grepl("entry withdrawn", dataRefs$Approved.name), ]
  withdrawnSymbols <- data.frame(Old=character(), New=character(), idx = numeric())
  
  for (idx in grep("withdrawn", dataRefs$Approved.symbol)){
    toAdd <- data.frame(str_extract(dataRefs$Approved.symbol[idx], ".*[^~withdrawn]"), word(dataRefs$Approved.name[idx], -1), idx)
    withdrawnSymbols <- rbind(withdrawnSymbols, toAdd)
    gsub("~withdrawn", "", dataRefs$Approved.symbol[idx])
  }
  rownames(dataRefs) <- dataRefs$Approved.symbol
  
  results <- list("withdrawnSymbols" = withdrawnSymbols, "HGNCData" = dataRefs)
  return(results)
}

# ====  PROCESS  ===============================================================
# Enter the step-by-step process of your project here. Strive to write your
# code so that you can simply run this entire block and re-create all
# intermediate results.
if (FALSE) {
  
  # ...
  
  
  
}

# ====  TESTS  =================================================================
if (FALSE) {
  # Enter your function tests here...
  
}


# [END]

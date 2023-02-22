## Create function that installs and loads the required packages 
## as specified in character vector (required.packages)


ipak <- function(required.packages){
  required.packages <- c(required.packages,"BiocManager")
  missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
  if(length(missing.packages)) {install.packages(missing.packages, dependencies = TRUE)}
  invisible(sapply(required.packages, library, character.only = TRUE))
}

ibiopak <- function(required.bioconductor.packages){
  #source("http://bioconductor.org/biocLite.R") ##replaced with "BiocManager"
  missing.packages <- required.bioconductor.packages[!(required.bioconductor.packages %in% installed.packages()[,"Package"])]
  #if(length(missing.packages)) {biocLite(missing.packages, dependencies = TRUE)}.
  if(length(missing.packages)) {BiocManager::install(missing.packages)}
  invisible(sapply(required.bioconductor.packages, library, character.only = TRUE))
}
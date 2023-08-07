rm(list=ls())

setwd("../")                                                                    # Setwd to location were output files have to be saved.

source("scripts/ipak.function.R")
source("scripts/get_functions.R")
library(data.table)
library(doRNG)

filename <- "BLCA"

outputpath <- paste0("Results/",filename)
sample_name <- paste0(filename,"_Full_");

#Get regulatory network
#=======================================================================================
rgbm_net <- read.table(gzfile(paste0(outputpath,"/Adjacency_Matrix/",sample_name,"Final_Adjacency_Matrix.csv.gz")),header=TRUE,sep=" ")

#Get the gene names and associate it with appropriate variables
#====================================================================================
load("Data/Others/TF_Target_Info.Rdata")
tfs <- tf_gene_names;
targets <- target_gene_names;
rownames(rgbm_net) <- tfs;
colnames(rgbm_net) <- targets;

#Make average regulon size to less than 100
#====================================================================================
while((sum(rgbm_net>0)/length(tfs))>100)
{
  median_weight <- median(rgbm_net[rgbm_net>0]);
  rgbm_net[rgbm_net<median_weight] <- 0;
}

#Get the edgelist for GENIE3
#####################################################################################
library(GENIE3)
set.seed(123)

#Get the expression matrix
z <- readRDS(file=paste0("Data/",filename,"/TCGA-",filename,"_normcounts.rda"))
rownames(z) <- targets
D <- as.matrix(log2(z+1))
rownames(D) <- rownames(z)
colnames(D) <- colnames(z)

#Build the genie network
genie_net <- GENIE3(D, nCores=12, targets=targets, regulators=tfs, returnMatrix=FALSE)
save(genie_net,file=paste0("Data/",filename,"/",filename,"_Genie.Rdata"))


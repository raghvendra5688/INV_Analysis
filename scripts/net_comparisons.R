rm(list=ls())

#setwd("../")                                                                    # Setwd to location were output files have to be saved.

source("scripts/ipak.function.R")
source("scripts/get_functions.R")
library(data.table)
library(doRNG)
library(igraph)
library(Matrix)

filenames <- c("BLCA","COAD","LGG","LUAD","KIRC","KIRP","MESO","PAAD","STAD","OV")

# for (filename in filenames)
# {
#   outputpath <- paste0("Results/",filename)
#   sample_name <- paste0(filename,"_Full_");
#   
#   #Get regulatory network
#   #=======================================================================================
#   rgbm_net <- read.table(gzfile(paste0(outputpath,"/Adjacency_Matrix/",sample_name,"Final_Adjacency_Matrix.csv.gz")),header=TRUE,sep=" ")
#   
#   #Get the gene names and associate it with appropriate variables
#   #====================================================================================
#   load("Data/Others/TF_Target_Info.Rdata")
#   tfs <- tf_gene_names;
#   targets <- target_gene_names;
#   rownames(rgbm_net) <- tfs;
#   colnames(rgbm_net) <- targets;
#   
#   #Make average regulon size to less than 100
#   #====================================================================================
#   while((sum(rgbm_net>0)/length(tfs))>100)
#   {
#     median_weight <- median(rgbm_net[rgbm_net>0]);
#     rgbm_net[rgbm_net<median_weight] <- 0;
#   }
#   
#   # #Get the edgelist for GENIE3
#   # #####################################################################################
#   # library(GENIE3)
#   # set.seed(123)
#   # 
#   # ##Get the expression matrix
#   # z <- readRDS(file=paste0("Data/",filename,"/TCGA-",filename,"_normcounts.rda"))
#   # rownames(z) <- targets
#   # D <- as.matrix(log2(z+1))
#   # rownames(D) <- rownames(z)
#   # colnames(D) <- colnames(z)
#   # 
#   # #Build the genie network
#   # genie_net <- GENIE3(D, nCores=36, targets=targets, regulators=tfs, returnMatrix=FALSE)
#   # save(genie_net,file=paste0("Data/",filename,"/",filename,"_Genie.Rdata"))
#   # 
#   
#   #Perform analysis on the genie networks
#   ##############################################################################
#   print(paste0("Currently at cancer: ",filename))
#   load(paste0("Data/",filename,"/",filename,"_Genie.Rdata"))
#   genie_network <- do.call(rbind,lapply(genie_net,matrix,ncol=length(tfs),byrow=TRUE))
#   genie_network[is.na(genie_network)]<- 0
#   rownames(genie_network) <- names(genie_net)
#   colnames(genie_network) <- names(genie_net[[1]])
#   
#   #Targets which are not in the GENIE network
#   missing_targets <- setdiff(targets,rownames(genie_network))
#   if (length(missing_targets)>0)
#   {
#     missing_target_matrix <- Matrix(0, nrow=length(missing_targets),ncol=length(colnames(genie_network)))
#     rownames(missing_target_matrix) <- missing_targets
#     colnames(missing_target_matrix) <- colnames(genie_network)
#     genie_network <- as.matrix(rbind(genie_network,missing_target_matrix))
#     colnames(genie_network) <- colnames(missing_target_matrix)
#     rownames(genie_network) <- c(names(genie_net),rownames(missing_target_matrix))
#   }
#   
#   #Make average regulon size to less than 100
#   #====================================================================================
#   while((sum(genie_network>0)/length(tfs))>100)
#   {
#     median_weight <- median(genie_network[genie_network>0]);
#     genie_network[genie_network<median_weight] <- 0;
#   }
#   
#   #Orient both RGBM and GENIE networks
#   ##############################################################################
#   rev_rgbm_net <- rgbm_net[tfs,targets]
#   rev_genie_network <- t(genie_network)[tfs,targets]
#   
#   save(rev_rgbm_net,file=paste0("Data/",filename,"/",filename,"_RGBM_Net.Rdata"))
#   save(rev_genie_network,file=paste0("Data/",filename,"/",filename,"_Genie_Net.Rdata"))  
#   rm(rgbm_net)
#   rm(genie_net)
#   rm(genie_network)
#   gc()
# }

#Once the small RGBM and GENIE networks are built, see similarity between them
###############################################################################
rgbm_no_edges <- list()
rgbm_specific_edges <- list()
genie_no_edges <- list()
genie_specific_edges <- list()
intersect_edges <- list()
union_edges <- list()
jc <- list()
for (i in c(1:length(filenames)))
{
  filename <- filenames[i]
  load(paste0("Data/",filename,"/",filename,"_RGBM_Net.Rdata"))
  load(paste0("Data/",filename,"/",filename,"_Genie_Net.Rdata"))
 
  rgbm_non_zero_ids <- which(rev_rgbm_net>0)
  genie_non_zero_ids <- which(rev_genie_network>0)
  
  intersect_edges[[i]] <- intersect(rgbm_non_zero_ids,genie_non_zero_ids)
  union_edges[[i]] <- union(rgbm_non_zero_ids,genie_non_zero_ids)
  jc[[i]] <- length(intersect_edges[[i]])/length(union_edges[[i]])
  
  rgbm_no_edges[[i]] <- length(rgbm_non_zero_ids)
  rgbm_specific_edges[[i]] <- length(setdiff(rgbm_non_zero_ids,genie_non_zero_ids))
  genie_no_edges[[i]] <- length(genie_non_zero_ids)
  genie_specific_edges[[i]] <- length(setdiff(genie_non_zero_ids,rgbm_non_zero_ids))
}

names(rgbm_no_edges) <- filenames
names(genie_no_edges) <- filenames
names(rgbm_specific_edges) <- filenames
names(genie_specific_edges) <- filenames
names(intersect_edges) <- filenames
names(union_edges) <- filenames
names(jc) <- filenames
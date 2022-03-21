library(data.table)
library(ggplot2)
library(doMC)
library(heatmap3)
library(gplots)
library(igraph)
library(yaGST)
library(parallel)
library(Matrix)
library(mltools)
require(viper)
require(limma)
require(GSVA)
registerDoMC(20)

setwd('../scripts/')

source('gene-reverse-network.R')
source('get_functions.R')

#Load mechanistic network
#======================================================================================
load('../Data/Others/me_net_full.Rdata')

#Load the RNASeq data
#=======================================================================================
filename <- "MESO"
outputpath <- paste0("../Results/ARACNE/",filename)
sample_name <- paste0(filename,"_Full_");
out <- loading_data(filename,M)
D <- as.matrix(log2(t(out[[1]])+1))

#Get regulatory network of ARACNE
#=======================================================================================
net <- read.table(paste0(outputpath,"/Adjacency_Matrix/Final_Adjacency_Matrix_v2.csv"),header=TRUE,sep=" ")

#Get the gene names and associate it with appropriate variables
#====================================================================================
tfs <- rownames(net)
targets <- colnames(net)

#Make average regulon size to less than 100
#====================================================================================
while((sum(net>0)/length(tfs))>100)
{
  median_weight <- median(net[net>0]);
  net[net<median_weight] <- 0;
}

#Convert the network into an edgelist for Viper
#====================================================================================
net_df <- get_viper_network(net, minsize=10)
write.table(net_df,paste0(outputpath,"/Adjacency_Matrix/Final_Adjacency_Matrix_for_Viper.csv"),row.names=F,col.names=F,quote=F,sep = "\t")


#Create the regulons for downstream enrichment analysis for Viper
#====================================================================================
afile <- paste0(outputpath,"/Adjacency_Matrix/Final_Adjacency_Matrix_for_Viper.csv")
viper_regulons <- aracne2regulon(afile, D, format=c("3col"))


#Get phenotype only for samples which have phenotype information
#=====================================================================================
load(paste0("../Data/",filename,"/",filename,"_INV_cluster_assignment_k2-6_v2.Rdata"))
high_low_output <- get_high_low_indices(table_cluster_assignment,D)
table_cluster_assignment <- high_low_output[[1]]
high_indices <- high_low_output[[2]]
low_indices <- high_low_output[[3]]
D <- high_low_output[[4]]
medium_indices <- setdiff(setdiff(c(1:dim(D)[2]),high_indices),low_indices)


#Get esets for ICR High and ICR Low and subtract baseline information from ICR Medium for TF-activity using Viper analysis
#=====================================================================================
D_high <- D[,high_indices]
D_low <- D[,low_indices]
D_medium <- D[,medium_indices]
vpsig_high <- viperSignature(D_high,D_medium,method="zscore",per=1000,seed=123,cores=20,verbose=TRUE)
vpres_high <- viper(vpsig_high, viper_regulons, minsize=10, cores=20, verbose = TRUE)
vpsig_low <- viperSignature(D_low,D_medium,method="zscore",per=1000,seed=123,cores=20,verbose=TRUE)
vpres_low <- viper(vpsig_low, viper_regulons, minsize=10, cores=20, verbose = TRUE)

vpres_all <- as.matrix(cbind(vpres_high,vpres_low))

#Perform limma to get the most differential gene sets
#=======================================================================================
inv_info <- c(rep("INV_High",length(high_indices)),rep("INV_Low",length(low_indices)))
design <- model.matrix(~ factor(inv_info))
fit <- lmFit(vpres_all, design)
fit <- eBayes(fit)
DEgeneSets <- topTable(fit, coef=2, number=Inf,
                       p.value=0.05, adjust="fdr")
save(DEgeneSets,file=paste0(outputpath,"/Adjacency_Matrix/",filename,"_All_Diff_MR_Info_ARACNE_Viper_BC_NES_1.Rdata"))

#Create the normalized enrichment matrix for the intersection genesets and plot it
#=================================================================================================
new_mr_activity_matrix <- as.matrix(vpres_all[rownames(DEgeneSets),])
rownames(new_mr_activity_matrix)
colnames(new_mr_activity_matrix) <- NULL

colcol <- rep("white",ncol(new_mr_activity_matrix));
high_ids <- c(1:length(high_indices));
low_ids <- c((length(high_indices)+1):length(colcol))
colcol[high_ids] <- "red"
colcol[low_ids] <- "blue"
hc_high.cols <- hclust(dist(t(new_mr_activity_matrix[,high_ids])),method='ward.D')
hc_low.cols <- hclust(dist(t(new_mr_activity_matrix[,low_ids])),method='ward.D')
hc.cols <- c(hc_high.cols$order,length(high_ids)+hc_low.cols$order);
new_mr_activity_matrix <- new_mr_activity_matrix[,hc.cols];

#Make the plot of normalized activity matrix clustered by hot vs cold immune response
#=================================================================================================
pdf(paste0(outputpath,"/",filename,"_Activity_ARACNE_Viper_BC_NES_1.pdf"),pointsize=9,height=10,width=10)
new_mr_activity_matrix[new_mr_activity_matrix <= quantile(new_mr_activity_matrix, 0.05)] <- quantile(new_mr_activity_matrix, 0.05)
new_mr_activity_matrix[new_mr_activity_matrix >= quantile(new_mr_activity_matrix, 0.95)] <- quantile(new_mr_activity_matrix, 0.95)
heatmap.2(x=new_mr_activity_matrix, Rowv = NULL, Colv = NULL, col = greenred(50), dendrogram = "none", trace="none", ColSideColors = colcol)
dev.off()
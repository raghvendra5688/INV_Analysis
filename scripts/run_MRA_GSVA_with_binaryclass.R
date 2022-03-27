library(data.table)
library(ggplot2)
library(doMC)
library(heatmap3)
library(gplots)
library(igraph)
library(yaGST)
library(parallel)
require(limma)
require(GSVA)
registerDoMC(20)

setwd('../scripts/')

source('gene-reverse-network.R')
source('get_functions.R')

#Load mechanistic network
#======================================================================================
load(gzfile('../Data/Others/me_net_full.Rdata.gz'))

#Load the RNASeq data
#=======================================================================================
filename <- "MESO"
outputpath <- paste0("../Results/",filename)
sample_name <- paste0(filename,"_Full_");
out <- loading_data(filename,M)
D <- as.matrix(log2(t(out[[1]])+1))

#Get regulatory network
#=======================================================================================
net <- read.table(gzfile(paste0(outputpath,"/Adjacency_Matrix/",sample_name,"Final_Adjacency_Matrix.csv.gz")),header=TRUE,sep=" ")

#Get the gene names and associate it with appropriate variables
#====================================================================================
load("../Data/Others/TF_Target_Info.Rdata")
tfs <- tf_gene_names;
targets <- target_gene_names;
rownames(net) <- tfs;
colnames(net) <- targets;

#Make average regulon size to less than 100
#====================================================================================
while((sum(net>0)/length(tfs))>100)
{
  median_weight <- median(net[net>0]);
  net[net<median_weight] <- 0;
}

#Get regulons for all TFs
#=====================================================================================
load(paste0(outputpath,"/Adjacency_Matrix/",sample_name,"regulon_FGSEA.Rdata"))

#Get phenotype only for samples which have phenotype information
#=====================================================================================
load(paste0("../Data/",filename,"/",filename,"_INV_cluster_assignment_k2-6_v2.Rdata"))
high_low_output <- get_high_low_indices(table_cluster_assignment,D)
table_cluster_assignment <- high_low_output[[1]]
high_indices <- high_low_output[[2]]
low_indices <- high_low_output[[3]]
D <- high_low_output[[4]]

#Perform gene set enrichment
#========================================================================================
gsvaout <- gsva(D, regulon_sets, verbose=FALSE, min.sz=10, max.sz = max(unlist(lapply(regulon_sets,length))), 
                mx.diff=TRUE, parallel.sz=20)

gsvaout_inv_related <- gsvaout[,c(high_indices,low_indices)]

#Perform limma to get the most differential gene sets
#=======================================================================================
inv_info <- c(rep("INV_High",length(high_indices)),rep("INV_Low",length(low_indices)))
design <- model.matrix(~ factor(inv_info))
fit <- lmFit(gsvaout_inv_related, design)
fit <- eBayes(fit)
DEgeneSets <- topTable(fit, coef=2, number=Inf,
                       p.value=0.05, adjust="BH")

FC <- 2^DEgeneSets$logFC
DEgeneSets$FC <- FC
save(DEgeneSets,file=paste0(outputpath,"/Adjacency_Matrix/",sample_name,"All_Diff_MR_Info_GSVA_BC_NES_1.Rdata"))

#Create the normalized enrichment matrix for the intersection genesets and plot it
#=================================================================================================
new_mr_activity_matrix <- as.matrix(gsvaout_inv_related[rownames(DEgeneSets),])
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
pdf(paste0(outputpath,"/Images/",sample_name,"Activity_GSVA_BC_NES_1.pdf"),pointsize=9,height=10,width=10)
new_mr_activity_matrix[new_mr_activity_matrix <= quantile(new_mr_activity_matrix, 0.05)] <- quantile(new_mr_activity_matrix, 0.05)
new_mr_activity_matrix[new_mr_activity_matrix >= quantile(new_mr_activity_matrix, 0.95)] <- quantile(new_mr_activity_matrix, 0.95)
heatmap.2(x=new_mr_activity_matrix, Rowv = NULL, Colv = NULL, col = greenred(50), dendrogram = "none", trace="none", ColSideColors = colcol)
dev.off()
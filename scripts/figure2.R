library(data.table)
library(ggplot2)
library(doMC)
library(heatmap3)
library(gplots)
library(igraph)
library(lattice)
library(parallel)
library(reshape2)
#source("https://bioconductor.org/biocLite.R")
#biocLite("fgsea")
require(fgsea)
library(fastcluster)
library(factoextra)
library(devtools)
library(imager)
library(gridExtra)
library(purrr)
library(dplyr)
library(jpeg)
library(png)
library(grid)
library(gridExtra)
library(grImport2)
library(ggpubr)
library(GEOquery)

#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

registerDoMC(20)

setwd('./scripts/')

source('gene-reverse-network.R')
source('get_functions.R')

#Load mechanistic network
#======================================================================================
load(gzfile(description = '../Data/Others/me_net_full.Rdata.gz'))

#Load the RNASeq data for BLCA
#=======================================================================================
filename <- "BLCA"
outputpath <- paste0("../Results/",filename)
sample_name <- paste0(filename,"_Full_");
out <- loading_data(filename,M)
D <- as.matrix(log2(t(out[[1]])+1))

#Filter to get the samples with the INV phenotype information
#======================================================================================
load(paste0("../Data/",filename,"/",filename,"_INV_cluster_assignment_k2-6_v2.Rdata"))
tcga_sample_ids <- rownames(table_cluster_assignment)
unmatched_tcga_sample_ids <- colnames(D);
filter_samples <- which(unmatched_tcga_sample_ids %in% tcga_sample_ids);
D <- D[,filter_samples];
jpeg(filename="../Results/Paper_Figures/Svgs_Jpgs/Supp_Figure_2.jpeg",width=900, height=480, units="px", pointsize=11)
boxplot(D[,1:100],ylab="Quantile-Normalized Counts",main="Quantile-Normalized Primary Tumor Samples")
dev.off()

order_indices <- NULL
for (i in 1:length(unmatched_tcga_sample_ids))
{
  order_indices <- c(order_indices,which(tcga_sample_ids==unmatched_tcga_sample_ids[i]));
}
table_cluster_assignment <- table_cluster_assignment[order_indices,];
high_indices <- which(table_cluster_assignment$HML_cluster=="INV High")
low_indices <- which(table_cluster_assignment$HML_cluster=="INV Low")
D1 <- D[,c(high_indices,low_indices)]
x <- apply(D1,1,IQR)
D2 <- D1[x>quantile(x,0.95),]

##Perform hierarchical clustering to get clusters on both rows and columns
hc.cols <- hclust(dist(t(D2)),method="ward.D2")
hc.rows <- hclust(dist(D2),method="ward.D2")

pdf(paste0("../Results/Paper_Figures/Figure2/RNA_Seq_Figure_2A.pdf"),width=6,height=8,pointsize = 11)
par(bg="white")
par(fg="black",col.axis="black",col.main="black",col.lab="white", cex.main=2.0)
heatmap.3(D2[hc.rows$order,hc.cols$order], Rowv = NULL, Colv = NULL, col = bluered(50), scale="none", main= "RNA-Seq Expression Data",
          xlab = "TCGA Samples", ylab = "Genes ", labRow = FALSE, labCol = FALSE, dendrogram = "none", cexRow = 6, cexCol = 6,
          key = TRUE, density.info = "none", KeyValueName = "log2(Normalized Expression)",
          margins = c(2,2), useRaster = TRUE)
dev.off()

###################################################################################################################

#Get regulatory network
#=======================================================================================
net <- read.table(gzfile(paste0(outputpath,"/Adjacency_Matrix/",sample_name,"Final_Adjacency_Matrix.csv.gz")),header=TRUE,sep=" ")

#Get the gene names and associate it with appropriate variables
#====================================================================================
load("../Data/Others/TF_Target_Info.Rdata")
tfs <- tf_gene_names;
targets <- target_gene_names;
rownames(net) <- tfs;
colnames(net) <- targets

#Get the network from cancer immunophenotype analysis
###################################################################################################################
#Load the TF-target signed regulon
#===================================================================================
#Get from cancer immunophenotype analysis

#######################################################################################################################
#Load the activity matrix for all TFs with regulon size >= 10
#==================================================================================
load(paste0(outputpath,"/Adjacency_Matrix/",sample_name,"Activity_matrix_FGSEA.Rdata"))
amat.pheno <- amat[,c(high_indices,low_indices)]
hc.cols <- hclust(dist(t(amat.pheno)),method="ward.D2")
hc.rows <- hclust(dist(amat.pheno),method="ward.D2")

min_activity <- min(amat.pheno)
max_activity <- max(amat.pheno)
for (i in 1:ncol(amat.pheno))
{
  for (j in 1:nrow(amat.pheno))
  {
    if (amat.pheno[j,i]>0)
    {
      amat.pheno[j,i] <- amat.pheno[j,i]/max_activity
    } else{
      amat.pheno[j,i] <- amat.pheno[j,i]/abs(min_activity)
    }
  }
  
}
hc.cols <- hclust(dist(t(amat.pheno)),method="average")
hc.rows <- hclust(dist(amat.pheno),method="ward.D2")

pdf("../Results/Paper_Figures/Figure2/BLCA_TF_Activity_Figure_2D.pdf",width=8,height=8,pointsize=11)
par(bg="white")
par(fg="black",col.axis="black",col.main="black",col.lab="black", cex.main=2.0, cex = 2.5)
heatmap.3(amat.pheno[hc.rows$order,hc.cols$order], Rowv = NULL, Colv = NULL, col = bluered(50), scale="none", main= "TF Regulon Activity Matrix",
          xlab = "TCGA Samples", ylab = "TFs", labRow = FALSE, labCol = FALSE, dendrogram = "none", cexRow=6, cexCol=6,
          key = TRUE, density.info = "none", KeyValueName = "Activity Value",
          margins = c(2,2), useRaster = TRUE)
dev.off()


###########################################################################################################################

#Get difference in expression of target genes using INV High vs INV Low phenotype
#===================================================================================
pheno_t <- c()
pv_t <- c()
for (j in 1:nrow(D)){
  cat(j,"\n")
  high_sens <- mean(D[j,high_indices]);
  low_sens <- mean(D[j,low_indices]);
  pheno_t <- c(pheno_t,high_sens-low_sens);
  pv_t <- c(pv_t,wilcox.test(D[j,high_indices],D[j,low_indices],exact=F)$p.value)
}
padj_t = p.adjust(pv_t,method="fdr")
names(pheno_t) <- rownames(D);
sorted_pheno_t <- sort(pheno_t,decreasing = TRUE)

#Perform gene set enrichment 
#========================================================================================
require(fgsea)
load(paste0(outputpath,"/Adjacency_Matrix/",sample_name,"regulon_FGSEA.Rdata"))
fgseaRes <- fgseaMultilevel(regulon_sets,sorted_pheno_t,minSize=10,maxSize=max(unlist(lapply(regulon_sets,length))))
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=15), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=15), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

#Get plot of top MRs and select top MRs
#========================================================================================
pdf("../Results/Paper_Figures/Figure2/BLCA_GeneRanks_FGSEA_Figure_2E.pdf",pointsize=11, height = 8, width = 8)
par(bg="white")
par(fg="black",col.axis="black",col.main="black",col.lab="black", cex.main=2.0, cex = 2.5)
plotGseaTable(regulon_sets[topPathways], sorted_pheno_t, fgseaRes, gseaParam = 0.5)
dev.off()

######################################################################################################################
#Get Top MR make their activity matrix
load(paste0(outputpath,"/Adjacency_Matrix/",sample_name,"TopMR_Info_FGSEA_BC_NES_1.Rdata"))
new_mr_activity_matrix <- as.matrix(amat[topmr_info$pathway,c(high_indices,low_indices)])
rownames(new_mr_activity_matrix)
colnames(new_mr_activity_matrix) <- NULL
min_activity <- min(new_mr_activity_matrix)
max_activity <- max(new_mr_activity_matrix)
for (i in 1:ncol(new_mr_activity_matrix))
{
  for (j in 1:nrow(new_mr_activity_matrix))
  {
    if (new_mr_activity_matrix[j,i]>0)
    {
      new_mr_activity_matrix[j,i] <- new_mr_activity_matrix[j,i]/max_activity
    } else{
      new_mr_activity_matrix[j,i] <- new_mr_activity_matrix[j,i]/abs(min_activity)
    }
  }
  
}

colcol <- matrix(0,nrow=ncol(new_mr_activity_matrix),ncol=1)
#rep("white",ncol(new_mr_activity_matrix));
high_ids <- c(1:length(high_indices));
low_ids <- c((length(high_indices)+1):length(colcol))
colcol[high_ids,1] <- "red"
colcol[low_ids,1] <- "blue"
colnames(colcol) <- "INV Phenotype"
hc_high.cols <- hclust(dist(t(new_mr_activity_matrix[,high_ids])),method='average')
hc_low.cols <- hclust(dist(t(new_mr_activity_matrix[,low_ids])),method='ward.D')
hc.cols <- c(hc_high.cols$order,length(high_ids)+hc_low.cols$order);
hc.rows <- hclust(dist(new_mr_activity_matrix),method='ward.D')

#Make the plot of activity matrix clustered by ICR High vs ICR Low
#=================================================================================================
pdf("../Results/Paper_Figures/Figure2/BLCA_TopMR_Activity_Figure_2F.pdf",height=8,width=8,pointsize=11)
par(bg="white")
par(fg="black",col.axis="black",col.main="black",col.lab="black", cex.main=2.0, cex = 2.5)
heatmap.3(new_mr_activity_matrix[hc.rows$order,hc.cols], Rowv = NULL, Colv = NULL, col = bluered(50), scale="none", main= "Master Regulator Activity Matrix",
          xlab = "TCGA Samples", ylab = "Top MRs", labRow = FALSE, labCol = FALSE, dendrogram = "none",
          key = TRUE, density.info = "none", KeyValueName = "Activity Value", ColSideColors = colcol, ColSideColorsSize=2,cexCol=2,
          margins = c(2,2), useRaster = TRUE)
dev.off()


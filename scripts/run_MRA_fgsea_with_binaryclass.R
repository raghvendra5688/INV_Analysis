library(data.table)
library(ggplot2)
library(doMC)
library(heatmap3)
library(gplots)
library(igraph)
library(yaGST)
library(parallel)
require(fgsea)
require(limma)
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

#Create correlation matrix
#====================================================================================
#corr_matrix <- cross.cor(D[tfs,],D,met="spearman",ncore = 20)
#corr_matrix[net==0] <- 0;
#save(corr_matrix,file=paste0(outputpath,"/Adjacency_Matrix/",sample_name,"Correlation_matrix.Rdata"))
load(gzfile(paste0(outputpath,"/Adjacency_Matrix/",sample_name,"Correlation_matrix.Rdata.gz")))

#Get regulons for all TFs
#=====================================================================================
regulon_sets <- get_regulons(net,corr_matrix,minsize=10)
save(regulon_sets,file=paste0(outputpath,"/Adjacency_Matrix/",sample_name,"regulon_FGSEA.Rdata"))
#load(paste0(outputpath,"/Adjacency_Matrix/",sample_name,"regulon_FGSEA.Rdata"))

#Get phenotype only for samples which have phenotype information
#=====================================================================================
load(paste0("../Data/",filename,"/",filename,"_INV_cluster_assignment_k2-6_v2.Rdata"))
high_low_output <- get_high_low_indices(table_cluster_assignment,D)
table_cluster_assignment <- high_low_output[[1]]
high_indices <- high_low_output[[2]]
low_indices <- high_low_output[[3]]
D <- high_low_output[[4]]

#Create activity matrix
#=====================================================================================
amat <- activity_mc(mexp = D, cormat = corr_matrix, tflist = names(regulon_sets), tau = 0.0)
save(amat,file=paste0(outputpath,"/Adjacency_Matrix/",sample_name,"Activity_matrix_FGSEA.Rdata"))
#load(paste0(outputpath,"/Adjacency_Matrix/",sample_name,"Activity_matrix_FGSEA.Rdata"))

#Get difference in expression of target genes using hot vs cold phenotype
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
fgseaRes <- fgseaMultilevel(regulon_sets,sorted_pheno_t,minSize=10,maxSize=max(unlist(lapply(regulon_sets,length))))
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=20), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=20), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

#Make plot of top MRs and select top MRs
#========================================================================================
pdf(paste0(outputpath,"/Images/",sample_name,"GeneRanks_FGSEA_BC.pdf"),pointsize=12,height = 10, width = 10)
plotGseaTable(regulon_sets[topPathways], sorted_pheno_t, fgseaRes, gseaParam = 0.5)
dev.off()

mra=as.data.frame(fgseaRes[order(pval)])
mra$leadingEdge=as.character(mra$leadingEdge)

#Get only those MRs with |NES|>1.0
topmr = mra$pathway[mra$padj<0.05 & abs(mra$NES)>1.0]
topmr_info = mra[mra$padj<0.05 & abs(mra$NES)>1.0,]
save(topmr_info,file=paste0(outputpath,"/Adjacency_Matrix/",sample_name,"TopMR_Info_FGSEA_BC_NES_1.Rdata"))

#Build the ordered activity matrix for significant MRs
#==============================================================================================
topmr <- as.character(as.vector(topmr_info$pathway))
amat.mr <- amat[topmr,];
amat.mr.high <- amat.mr[,high_indices];
amat.mr.low <- amat.mr[,low_indices]
wilcox_test_info <- perform_wilcox_test(amat.mr.high,amat.mr.low)
wilcox_test_info <- wilcox_test_info[order(wilcox_test_info$FC_Mean,decreasing=T),]
write.table(wilcox_test_info,file=paste0(outputpath,"/Adjacency_Matrix/",sample_name,"WilCox_Test_FGSEA_MR_Result_NES_1.csv"),row.names=T,col.names=T);

#Prepare the activity matrix to plot
#=================================================================================================
new_mr_activity_matrix <- as.matrix(amat[rownames(wilcox_test_info),c(high_indices,low_indices)])
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

#Make the plot of activity matrix clustered by hot vs cold immune response
#=================================================================================================
pdf(paste0(outputpath,"/Images/",sample_name,"Activity_FGSEA_BC_NES_1.pdf"),pointsize=9,height=10,width=10)
new_mr_activity_matrix[new_mr_activity_matrix <= quantile(new_mr_activity_matrix, 0.05)] <- quantile(new_mr_activity_matrix, 0.05)
new_mr_activity_matrix[new_mr_activity_matrix >= quantile(new_mr_activity_matrix, 0.95)] <- quantile(new_mr_activity_matrix, 0.95)
heatmap.2(x=new_mr_activity_matrix, Rowv = NULL, Colv = NULL, col = greenred(50), dendrogram = "none", trace="none", ColSideColors = colcol)
dev.off()
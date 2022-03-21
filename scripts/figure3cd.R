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
library(plyr)
library(RColorBrewer)
library(ggsignif)
library(gridExtra)
library(ggpubr)
library(grImport2)
library(matrixStats)
warnings("off")

source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

registerDoMC(20)

setwd('../scripts/')

source('gene-reverse-network.R')
source('get_functions.R')

revise_wilcox_test_info <- function(wilcox_test_info)
{
  wilcox_test_info$Pval <- round(wilcox_test_info$Pval,5)
  wilcox_test_info$Mean1 <- round(wilcox_test_info$Mean1,3)
  wilcox_test_info$Mean2 <- round(wilcox_test_info$Mean2,3)
  wilcox_test_info$FC_Mean <- round(wilcox_test_info$FC_Mean,3)
  return(wilcox_test_info)
}

all_inv <- c("LGG","KIRP","PAAD","MESO","KIRC",
             "COAD","BLCA","STAD","LUAD","OV")

common_mrs_all_list <- NULL

output_df <- read.table("../Results/Paper_Text/All_TF_Activity_Information_v2.csv",header = TRUE)
output_df$MR <- as.character(as.vector(output_df$MR))
output_df$Cancer <- as.character(as.vector(output_df$Cancer))
colors <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#a9a9a9', '#008080', '#e6beff', '#00000f','#f00000')

common_mrs_all_list <- c(output_df[output_df$Cancer==all_inv[1],]$MR)
for (i in 2:length(all_inv))
{
  new_mrs <- output_df[output_df$Cancer==all_inv[i],]$MR
  common_mrs_all_list <- intersect(common_mrs_all_list,new_mrs)
}

common_mrs_all_activity_matrix <- matrix(0,nrow=length(common_mrs_all_list),ncol=2*length(all_inv))
common_mrs_all_list <- sort(common_mrs_all_list)
rownames(common_mrs_all_activity_matrix) <- common_mrs_all_list
colnames(common_mrs_all_activity_matrix) <- c(all_inv, all_inv)
for (i in 1:length(all_inv))
{
  cancer <- all_inv[i]
  temp_df <- output_df[output_df$MR %in% common_mrs_all_list & output_df$Cancer==cancer,]
  temp_df <- temp_df[order(temp_df$MR),]
  inv_high <- temp_df$Median_INV_High
  inv_low <- temp_df$Median_INV_Low
  common_mrs_all_activity_matrix[common_mrs_all_list,i] <- inv_high
  common_mrs_all_activity_matrix[common_mrs_all_list,(i+length(all_inv))] <- inv_low
}

colcol <- matrix(0,nrow=2*length(all_inv),ncol=1)
colcol[c(1:length(all_inv)),1] <- "red"
colcol[c((length(all_inv)+1):(2*length(all_inv))),1] <- "blue"

#Figure 3C 
#=========================================================================================================
pdf("../Results/Paper_Figures/Figure3/Common_TopMRs_Activity_INV_Figure_3C.pdf",height = 10, width=12, pointsize = 14)
heatmap.3(common_mrs_all_activity_matrix, Rowv = TRUE, Colv=TRUE, col = bluered(100), scale="none", main= "Median Activity of MRs common across all Invasive Cancers",
          labRow = rownames(common_mrs_all_activity_matrix), labCol = colnames(common_mrs_all_activity_matrix), dendrogram = "col", 
          key = TRUE, density.info = "none", KeyValueName = "Activity Value", ColSideColors = colcol, ColSideColorsSize = 2,
          margins = c(6,6), useRaster = FALSE, cexRow = 1.25, cexCol = 2, cellnote = ifelse(common_mrs_all_activity_matrix>0,"+","-"), notecex = 2, notecol = "black")
dev.off()

#Supplementary information about top common MRs for INVasive cancer based on fold change values
A <- common_mrs_all_activity_matrix[common_mrs_all_list,c(1:length(all_inv))];
B <- common_mrs_all_activity_matrix[common_mrs_all_list,c((length(all_inv)+1):(2*length(all_inv)))]
wilcox_test_info <- perform_wilcox_test(A,B)
wilcox_test_info <- wilcox_test_info[order(wilcox_test_info$FC_Mean,decreasing=TRUE),]
wilcox_test_info <- revise_wilcox_test_info(wilcox_test_info)
write.table(wilcox_test_info,"../Results/Paper_Text/Supplementary_Table_S1_Common_MRs_INV.csv",row.names=T,col.names=T,quote=F)

inv_activity_df <- NULL
for (cancer in all_inv)
{
  load(paste0("../Results/",cancer,"/Adjacency_Matrix/",cancer,"_Full_Activity_matrix_FGSEA.Rdata"))
  amat[amat>0] <- amat[amat>0]/max(amat)
  amat[amat<0] <- amat[amat<0]/abs(min(amat))
  
  #Load mechanistic network
  #======================================================================================
  load('../Data/Others/me_net_full.Rdata')
  
  #Get the data
  filename = cancer
  out <- loading_data(filename,M)
  D <- as.matrix(log2(t(out[[1]])+1))
  
  #Get high and low indices
  load(paste0("../Data/",filename,"/",filename,"_INV_cluster_assignment_k2-6_v2.Rdata"))
  high_low_output <- get_high_low_indices(table_cluster_assignment,D)
  table_cluster_assignment <- high_low_output[[1]]
  high_indices <- high_low_output[[2]]
  low_indices <- high_low_output[[3]]
  
  for (topmr in common_mrs_all_list)
  {
    high_activity <- amat[topmr,high_indices]
    low_activity <- amat[topmr,low_indices]
    inv_activity_df <- rbind(inv_activity_df,
                                 cbind(rep(cancer,length(high_activity)+length(low_activity)),
                                       rep(topmr,length(high_activity)+length(low_activity)),
                                       c(rep("INV High",length(high_activity)),rep("INV Low",length(low_activity))),
                                       c(high_activity,low_activity)))
  }
}

inv_activity_df <- as.data.frame(inv_activity_df)
colnames(inv_activity_df) <- c("Cancer","TopMR","Phenotype","Activity_Value")
inv_activity_df$Cancer <- as.character(as.vector(inv_activity_df$Cancer))
inv_activity_df$TopMR <- as.character(as.vector(inv_activity_df$TopMR))
inv_activity_df$Phenotype <- as.character(as.vector(inv_activity_df$Phenotype))
inv_activity_df$Activity_Value <- as.numeric(as.vector(inv_activity_df$Activity_Value))

p3 <- ggplot(data = inv_activity_df, aes(x=TopMR, y=Activity_Value)) + 
  geom_boxplot(aes(fill=Phenotype)) + facet_wrap( ~ Cancer, nrow=2, ncol=5) + xlab("Top Common MRs") + ylab("Activity Value") +
  geom_point(aes(y=Activity_Value, group=Phenotype), size=0.01, position = position_dodge(width=0.75))+
  guides(fill=guide_legend(title="Phenotype")) + scale_fill_manual(values=c("red","blue")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.text.x = element_text(angle = -90)) +
  ggtitle("Activities of Common MRs common across Invasive Cancers ") + theme(text = element_text(size=12)) + theme(plot.title = element_text(hjust = 0.5))

ggsave(file="../Results/Paper_Figures/Figure3/Common_TopMRs_Activity_Figure_3D.pdf",plot = p3, device = pdf(), width = 12, height=12, units = "in", dpi = 300)
dev.off()
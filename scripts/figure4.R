library(data.table)
library(ggplot2)
library(doMC)
library(heatmap3)
library(gplots)
library(igraph)
library(lattice)
library(parallel)
library(reshape2)
require(fgsea)
library(ComplexHeatmap)
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
library(plotly)
library(processx)
library(ggpubr)
library(grImport2)
#library(clusterProfiler)
library(magrittr)
library(msigdbr)
library(networkD3)
library(matrixStats)
library(rbokeh)
warnings("off")

source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

registerDoMC(20)

setwd('../scripts/')

source('gene-reverse-network.R')
source('get_functions.R')


get_tf_median_activity <- function(cancer)
{
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
  
  load(paste0("../Results/",cancer,"/Adjacency_Matrix/",cancer,"_Full_Activity_matrix_FGSEA.Rdata"))
  amat[amat>0] <- amat[amat>0]/max(amat)
  amat[amat<0] <- amat[amat<0]/abs(min(amat))
  
  tfs_median_high_activity <- rowMedians(amat[,high_indices])
  tfs_median_low_activity <- rowMedians(amat[,low_indices])
  activity_info <- cbind(rownames(amat),rep(cancer,nrow(amat)),tfs_median_high_activity,tfs_median_low_activity)
  return(activity_info)
}

get_tf_activity_high_low <- function(cancer,tf)
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
  
  
  if (tf %in% rownames(amat))
  {
    tf_high_activity <- as.numeric(as.vector(amat[tf,high_indices]))
    tf_low_activity <- as.numeric(as.vector(amat[tf,low_indices]))
  }
  else{
    tf_high_activity <- rep(0,length(high_indices))
    tf_low_activity <- rep(0,length(low_indices))
  }
  return(list(tf_high_activity,tf_low_activity))
}


get_common_mrs <- function(cancer)
{
  load(paste0("../Results/",cancer,"/Adjacency_Matrix/",cancer,"_Full_TopMR_Info_FGSEA_BC_NES_1.Rdata"))
  load(paste0("../Results/",cancer,"/Combined/",cancer,"_Common_TopMRs.Rdata"))
  common_topmr <- topmr_info[topmr_info$pathway %in% common_mrs,]$pathway
  list_common_mrs <- as.character(as.vector(common_topmr));
  return(list_common_mrs)
}

INV_enabled <- c("BLCA","COAD","LGG","OV","PRAD","PAAD","STAD","LUAD","LUSC","ESCA")
INV_disabled <- c("HNSC","DLBC","LIHC","THCA")

tfs_enabled_activity_info <- NULL
tfs_disabled_activity_info <- NULL

#Get the list of TFs which are common to all the 12 ICR Enabled cancers and their median activity in ICR High vs ICR Low samples per cancer
#############################################################################################################################
colors <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#a9a9a9', '#008080', '#e6beff', '#00000f','#f00000')
first_cancer <- INV_enabled[1];
tfs_enabled_activity_info <- get_tf_median_activity(first_cancer)
for (i in 2:length(INV_enabled))
{
  cancer <- INV_enabled[i]
  temp <- get_tf_median_activity(cancer)
  tfs_enabled_activity_info <- rbind(tfs_enabled_activity_info,temp)
}
tfs_enabled_activity_info <- as.data.frame(tfs_enabled_activity_info)
colnames(tfs_enabled_activity_info) <- c("TFs","Cancer","Median_INV_High","Median_INV_Low")
tfs_enabled_activity_info$TFs <- as.character(as.vector(tfs_enabled_activity_info$TFs))
tfs_enabled_activity_info$Cancer <- as.character(as.vector(tfs_enabled_activity_info$Cancer))
tfs_enabled_activity_info$Median_INV_High <- as.numeric(as.vector(tfs_enabled_activity_info$Median_INV_High))
tfs_enabled_activity_info$Median_INV_Low <- as.numeric(as.vector(tfs_enabled_activity_info$Median_INV_Low))

#Select those transcription regulators which are TFs with regulons of size >= 10 in all the Inv Enabled Cancers
common_tfs_enabled_list <- names(which(table(tfs_enabled_activity_info$TFs)>9))
common_tfs_enabled_activity_info <- tfs_enabled_activity_info[tfs_enabled_activity_info$TFs %in% common_tfs_enabled_list,]

#Make a list of how many times does a common MR appear for the INV Enabled cancers
first_cancer <- INV_enabled[1];
list_common_mrs_enabled <- get_common_mrs(first_cancer)
topMRs_enabled <- rep(list("BLCA"),length(list_common_mrs_enabled))
names(topMRs_enabled) <- list_common_mrs_enabled
for(i in 2:length(INV_enabled))
{
  cancer_type <- INV_enabled[i]
  list_common_mrs_enabled <- get_common_mrs(cancer_type)
  for (j in 1:length(list_common_mrs_enabled))
  {
    if (list_common_mrs_enabled[j] %in% names(topMRs_enabled))
    {
      topMRs_enabled[[list_common_mrs_enabled[j]]] <- c(topMRs_enabled[[list_common_mrs_enabled[j]]],cancer_type)
    } else {
      topMRs_enabled[[list_common_mrs_enabled[j]]] <- cancer_type
    }
  }
}
all_icr_cancers <- unlist(lapply(topMRs_enabled,function(x) paste(x,collapse=" ")))
no_icr_cancers <- unlist(lapply(topMRs_enabled,function(x) length(x)))
enabled_df <- cbind(names(all_icr_cancers),as.character(all_icr_cancers),as.numeric(no_icr_cancers))
enabled_df <- as.data.frame(enabled_df)
colnames(enabled_df) <- c("MR","List_INV_Cancers","No_INV_Cancers")
enabled_df$MR <- as.character(as.vector(enabled_df$MR))
enabled_df$List_INV_Cancers <- as.character(as.vector(enabled_df$List_INV_Cancers))
enabled_df$No_INV_Cancers <- as.character(as.vector(enabled_df$No_INV_Cancers))
enabled_df <- enabled_df[order(enabled_df$No_INV_Cancers,decreasing=T),]

#Selct MRs which are present in atleast 50% of the INV enabled cancers and are also TF for all INV enabled cancers
#enabled_7_MR <- enabled_df[enabled_df$No_ICR_Cancers==7,]$MR
#enabled_7_MR_and_TF <- enabled_7_MR[enabled_7_MR %in% common_tfs_enabled_list]
#enabled_7_MR_and_TF <- paste0(enabled_7_MR_and_TF,collapse = ", ")
temp_enabled_df <- enabled_df[enabled_df$No_INV_Cancers>=5,]
temp_enabled_df <- temp_enabled_df[temp_enabled_df$MR %in% common_tfs_enabled_list, ]
write.table(temp_enabled_df,"../Results/Paper_Text/All_MRs_TFs_INV_Enabled_Supplementary_Table_S3.csv",row.names=F,col.names=T,quote=F,sep=",")

#Get list of common MRs which are all TFs in 5 out of 10 cancers
#########################################################################################################################
common_mrs_enabled_list <- temp_enabled_df$MR;
common_mrs_enabled_list <- sort(common_mrs_enabled_list)
common_mrs_enabled_activity_matrix <- matrix(0,nrow=length(common_mrs_enabled_list),ncol=2*length(INV_enabled))
rownames(common_mrs_enabled_activity_matrix) <- common_mrs_enabled_list
colnames(common_mrs_enabled_activity_matrix) <- c(INV_enabled,INV_enabled)
for (i in 1:length(INV_enabled))
{
  cancer <- INV_enabled[i]
  temp_df <- common_tfs_enabled_activity_info[common_tfs_enabled_activity_info$TFs %in% common_mrs_enabled_list 
                                              & common_tfs_enabled_activity_info$Cancer==cancer,]
  temp_df <- temp_df[order(temp_df$TFs),]
  inv_high <- temp_df$Median_INV_High
  inv_low <- temp_df$Median_INV_Low
  common_mrs_enabled_activity_matrix[common_mrs_enabled_list,i] <- inv_high
  common_mrs_enabled_activity_matrix[common_mrs_enabled_list,(i+length(INV_enabled))] <- inv_low
}

colcol <- matrix(0,nrow=2*length(INV_enabled),ncol=1)
colcol[c(1:length(INV_enabled)),1] <- "yellow"
colcol[c((length(INV_enabled)+1):(2*length(INV_enabled))),1] <- "green"

#Perform Wilcox ranksum test or Mann-Whitney test to identify the MRs whose activity between the INV High and INV Low are significant for INV Enabled cancers
#######################################################################################################
output_df <- common_tfs_enabled_activity_info
nes_info <- NULL
nes_pval <- NULL
for (i in 1:length(common_mrs_enabled_list))
{
  mr <- common_mrs_enabled_list[i]
  print(paste0("Stuck at MR: ",mr," which is MR no: ",i))
  mr_inv_high_values <- NULL
  mr_inv_low_values <- NULL
  for (cancer in INV_enabled)
  {
    mr_activity_list <- get_tf_activity_high_low(cancer,mr)
    mr_inv_high_values <- c(mr_inv_high_values,mr_activity_list[[1]])
    mr_inv_low_values <- c(mr_inv_low_values,mr_activity_list[[2]])
  }
  nes_info <- c(nes_info,median(mr_inv_high_values)-median(mr_inv_low_values))
  nes_pval <- c(nes_pval,wilcox.test(mr_inv_high_values,mr_inv_low_values,exact=F)$p.value)
}
inv_enabled_median_comparison <- cbind(common_mrs_enabled_list,nes_info,nes_pval)
inv_enabled_median_comparison <- as.data.frame(inv_enabled_median_comparison)
colnames(inv_enabled_median_comparison) <- c("MR","FC_Median","Pval")
inv_enabled_median_comparison$MR <- as.character(as.vector(inv_enabled_median_comparison$MR))
inv_enabled_median_comparison$FC_Median <- round(as.numeric(as.vector(inv_enabled_median_comparison$FC_Median)),3)
inv_enabled_median_comparison$Pval <- p.adjust(as.numeric(as.vector(inv_enabled_median_comparison$Pval)),method = "fdr")
final_common_mrs_enabled_list <- inv_enabled_median_comparison[inv_enabled_median_comparison$Pval<=0.05,]$MR

final_common_mrs_enabled_activity_matrix <- common_mrs_enabled_activity_matrix[final_common_mrs_enabled_list,]

#Figure 4A 
#=========================================================================================================
pdf("../Results/Paper_Figures/Figure4/Common_TopMRs_Activity_INV_Enabled_Figure_4A.pdf",height = 12, width=14, pointsize = 14)
par(bg="white")
par(fg="black",col.axis="black",col.main="black",col.lab="black", cex.main=1.75)
p1 <- heatmap.3(final_common_mrs_enabled_activity_matrix, Rowv = TRUE, Colv=, col = bluered(100), scale="none", main= "Median Activity of Common MRs in INV Enabled Cancers", # (>=4 out of 8)",
                dendrogram = "both", key = TRUE, density.info = "none", KeyValueName = "Activity Value", ColSideColors = colcol, ColSideColorsSize = 2,
                margins = c(6,6), useRaster = FALSE, cexRow = 0.6, cexCol = 2.0, cellnote = ifelse(final_common_mrs_enabled_activity_matrix>0,"+","-"), notecex = 1, notecol = "black")
dev.off()

#Identify the MRs which are specific to ICR Low
ordered_mrs_enabled <- rev(rownames(final_common_mrs_enabled_activity_matrix)[p1$rowInd])
ordered_p_values_inv_enabled <- NULL
for (i in 1:length(ordered_mrs_enabled))
{
  mr <- ordered_mrs_enabled[i]
  adj_pval <- inv_enabled_median_comparison[inv_enabled_median_comparison$MR==mr,]$Pval
  temp <- cbind(mr,adj_pval)
  ordered_p_values_inv_enabled <- rbind(ordered_p_values_inv_enabled,temp)
}
ordered_p_values_inv_enabled <- as.data.frame(ordered_p_values_inv_enabled)
colnames(ordered_p_values_inv_enabled) <- c("MR","Adj_Pval")
ordered_p_values_inv_enabled$MR <- as.character(as.vector(ordered_p_values_inv_enabled$MR))
ordered_p_values_inv_enabled$Adj_Pval <- as.numeric(as.vector(ordered_p_values_inv_enabled$Adj_Pval))
write.table(ordered_p_values_inv_enabled,"../Results/Paper_Text/Ordered_Pvalues_INV_Enabled_MRs.csv",row.names=F,col.names=T,quote=F,sep=",")

#Get activity of MRs specific to INV Low of interest from enabled cancers
mrs_of_interest_enabled <- inv_enabled_median_comparison[inv_enabled_median_comparison$FC_Median<0 &
                                                           inv_enabled_median_comparison$Pval<=0.05,]$MR
enabled_activity_df <- NULL
for (cancer in INV_enabled)
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
  
  
  for (topmr in mrs_of_interest_enabled)
  {
    high_activity <- amat[topmr,high_indices]
    low_activity <- amat[topmr,low_indices]
    enabled_activity_df <- rbind(enabled_activity_df,
                                 cbind(rep(cancer,length(high_activity)+length(low_activity)),
                                       rep(topmr,length(high_activity)+length(low_activity)),
                                       c(rep("INV High",length(high_activity)),rep("INV Low",length(low_activity))),
                                       c(high_activity,low_activity)))
  }
}

enabled_activity_df <- as.data.frame(enabled_activity_df)
colnames(enabled_activity_df) <- c("Cancer","TopMR","Phenotype","Activity_Value")
enabled_activity_df$Cancer <- as.character(as.vector(enabled_activity_df$Cancer))
enabled_activity_df$TopMR <- as.character(as.vector(enabled_activity_df$TopMR))
enabled_activity_df$Phenotype <- as.character(as.vector(enabled_activity_df$Phenotype))
enabled_activity_df$Activity_Value <- as.numeric(as.vector(enabled_activity_df$Activity_Value))

supp_p3 <- ggplot(data = enabled_activity_df, aes(x=Cancer, y=Activity_Value)) + 
  geom_boxplot(aes(fill=Phenotype)) + facet_wrap( ~ TopMR, nrow=11, ncol=10) + xlab("Cancer") + ylab("Activity Value") +
  geom_point(aes(y=Activity_Value, group=Phenotype), size=0.01, position = position_dodge(width=0.75))+
  guides(fill=guide_legend(title="Phenotype")) + scale_fill_manual(values=c("yellow","green")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.text.x = element_text(angle = -90)) +
  ggtitle("Activities of Top MRs of specific to INV Low in INV Enabled Cancers (>=5 out of 10)") + theme(text = element_text(size=12)) + theme(plot.title = element_text(hjust = 0.5))
ggsave(file="../Results/Paper_Figures/Svgs_Jpgs/Supp_Figure_4.jpg",plot = supp_p3, device = jpeg(), width = 14, height=12, units = "in", dpi = 300)
dev.off()

inv_enabled_median_comparison <- inv_enabled_median_comparison[order(inv_enabled_median_comparison$FC_Median,decreasing = T),]
final_inv_enabled_median_comparison <- inv_enabled_median_comparison[inv_enabled_median_comparison$Pval<=0.05,]
final_inv_enabled_median_comparison$Pval <- as.numeric(as.vector(signif(final_inv_enabled_median_comparison$Pval,digits=6)))
write.table(final_inv_enabled_median_comparison,file="../Results/Paper_Text/All_MRs_TFs_INV_Enabled_Supplementary_Table_S4.csv",
            row.names = F, col.names = T, sep = " & ", quote=F)


#Perform Analysis for INV Disabled
###########################################################################################################################
tfs_disabled_activity_info <- get_tf_median_activity(INV_disabled[1])
for (i in 2:length(INV_disabled))
{
  cancer <- INV_disabled[i]
  temp <- get_tf_median_activity(cancer)
  tfs_disabled_activity_info <- rbind(tfs_disabled_activity_info,temp)
}
tfs_disabled_activity_info <- as.data.frame(tfs_disabled_activity_info)
colnames(tfs_disabled_activity_info) <- c("TFs","Cancer","Median_INV_High","Median_INV_Low")
tfs_disabled_activity_info$TFs <- as.character(as.vector(tfs_disabled_activity_info$TFs))
tfs_disabled_activity_info$Cancer <- as.character(as.vector(tfs_disabled_activity_info$Cancer))
tfs_disabled_activity_info$Median_INV_High <- as.numeric(as.vector(tfs_disabled_activity_info$Median_INV_High))
tfs_disabled_activity_info$Median_INV_Low <- as.numeric(as.vector(tfs_disabled_activity_info$Median_INV_Low))

#Select those transcription regulators which are TFs with regulons of size >= 10 in all the INV Disabled Cancers
common_tfs_disabled_list <- names(which(table(tfs_disabled_activity_info$TFs)>3))
common_tfs_disabled_activity_info <- tfs_disabled_activity_info[tfs_disabled_activity_info$TFs %in% common_tfs_disabled_list,]

first_cancer <- INV_disabled[1];
list_common_mrs_disabled <- get_common_mrs(first_cancer)
topMRs_disabled <- rep(list("HNSC"),length(list_common_mrs_disabled))
names(topMRs_disabled) <- list_common_mrs_disabled
for(i in 2:length(INV_disabled))
{
  cancer_type <- INV_disabled[i]
  list_common_mrs_disabled <- get_common_mrs(cancer_type)
  for (j in 1:length(list_common_mrs_disabled))
  {
    if (list_common_mrs_disabled[j] %in% names(topMRs_disabled))
    {
      topMRs_disabled[[list_common_mrs_disabled[j]]] <- c(topMRs_disabled[[list_common_mrs_disabled[j]]],cancer_type)
    } else {
      topMRs_disabled[[list_common_mrs_disabled[j]]] <- cancer_type
    }
  }
}

all_inv_cancers <- unlist(lapply(topMRs_disabled,function(x) paste(x,collapse=" ")))
no_inv_cancers <- unlist(lapply(topMRs_disabled,function(x) length(x)))
disabled_df <- cbind(names(all_inv_cancers),as.character(all_inv_cancers),as.numeric(no_inv_cancers))
disabled_df <- as.data.frame(disabled_df)
colnames(disabled_df) <- c("MR","List_INV_Cancers","No_INV_Cancers")
disabled_df$MR <- as.character(as.vector(disabled_df$MR))
disabled_df$List_INV_Cancers <- as.character(as.vector(disabled_df$List_INV_Cancers))
disabled_df$No_INV_Cancers <- as.character(as.vector(disabled_df$No_INV_Cancers))
disabled_df <- disabled_df[order(disabled_df$No_INV_Cancers,decreasing=T),]

temp_disabled_df <- disabled_df[disabled_df$No_INV_Cancers>=2,]
temp_disabled_df <- temp_disabled_df[temp_disabled_df$MR %in% common_tfs_disabled_list, ]
write.table(temp_disabled_df,"../Results/Paper_Text/All_MRs_TFs_INV_Disabled_Supplementary_Table_S5.csv",row.names=F,col.names=T,quote=F,sep=",")

#Get list of common MRs which are also TFs in 2 out of 4 cancers
############################################################################################
common_mrs_disabled_list <- temp_disabled_df$MR;
common_mrs_disabled_list <- sort(common_mrs_disabled_list)
common_mrs_disabled_activity_matrix <- matrix(0,nrow=length(common_mrs_disabled_list),ncol=2*length(INV_disabled))
rownames(common_mrs_disabled_activity_matrix) <- common_mrs_disabled_list
colnames(common_mrs_disabled_activity_matrix) <- c(INV_disabled,INV_disabled)
for (i in 1:length(INV_disabled))
{
  cancer <- INV_disabled[i]
  temp_df <- common_tfs_disabled_activity_info[common_tfs_disabled_activity_info$TFs %in% common_mrs_disabled_list 
                                               & common_tfs_disabled_activity_info$Cancer==cancer,]
  temp_df <- temp_df[order(temp_df$TFs),]
  inv_high <- temp_df$Median_INV_High
  inv_low <- temp_df$Median_INV_Low
  common_mrs_disabled_activity_matrix[common_mrs_disabled_list,i] <- inv_high
  common_mrs_disabled_activity_matrix[common_mrs_disabled_list,(i+length(INV_disabled))] <- inv_low
}

colcol <- matrix(0,nrow=2*length(INV_disabled),ncol=1)
colcol[c(1:length(INV_disabled)),1] <- "yellow"
colcol[c((length(INV_disabled)+1):(2*length(INV_disabled))),1] <- "green"

#Perform Wilcox ranksum test or Mann-Whitney test to identify the MRs whose activity between the INV High and INV Low are significant for INV Disabled cancers
#######################################################################################################
output_df <- common_tfs_disabled_activity_info
nes_info <- NULL
nes_pval <- NULL
for (i in 1:length(common_mrs_disabled_list))
{
  mr <- common_mrs_disabled_list[i]
  mr_inv_high_values <- NULL
  mr_inv_low_values <- NULL
  for (cancer in INV_disabled)
  {
    mr_activity_list <- get_tf_activity_high_low(cancer,mr)
    mr_inv_high_values <- c(mr_inv_high_values,mr_activity_list[[1]])
    mr_inv_low_values <- c(mr_inv_low_values,mr_activity_list[[2]])
  }
  nes_info <- c(nes_info,median(mr_inv_high_values)-median(mr_inv_low_values))
  nes_pval <- c(nes_pval,wilcox.test(mr_inv_high_values,mr_inv_low_values,exact=F)$p.value)
}
inv_disabled_median_comparison <- cbind(common_mrs_disabled_list,nes_info,nes_pval)
inv_disabled_median_comparison <- as.data.frame(inv_disabled_median_comparison)
colnames(inv_disabled_median_comparison) <- c("MR","FC_Median","Pval")
inv_disabled_median_comparison$MR <- as.character(as.vector(inv_disabled_median_comparison$MR))
inv_disabled_median_comparison$FC_Median <- round(as.numeric(as.vector(inv_disabled_median_comparison$FC_Median)),3)
inv_disabled_median_comparison$Pval <- p.adjust(as.numeric(as.vector(inv_disabled_median_comparison$Pval)),method = "fdr")
final_common_mrs_disabled_list <- inv_disabled_median_comparison[inv_disabled_median_comparison$Pval<=0.05,]$MR
final_common_mrs_disabled_activity_matrix <- common_mrs_disabled_activity_matrix[final_common_mrs_disabled_list,]

#Figure 4B
#=========================================================================================================
pdf("../Results/Paper_Figures/Figure4/Common_TopMRs_Activity_INV_Disabled_Figure_4B.pdf",height = 12, width=14, pointsize = 14)
par(bg="white")
par(fg="black",col.axis="black",col.main="black",col.lab="black", cex.main=1.65)
p2 <- heatmap.3(final_common_mrs_disabled_activity_matrix, Rowv = TRUE, Colv=TRUE, col = bluered(100), scale="none", main= "Median Activity of Common MRs in INV Disabled Cancers", # (>=2 out of 4)",
                dendrogram = "both", key = TRUE, density.info = "none", KeyValueName = "Activity Value", ColSideColors = colcol, ColSideColorsSize = 2,
                margins = c(6,6), useRaster = FALSE, cexRow = 0.6, cexCol = 2, cellnote = ifelse(final_common_mrs_disabled_activity_matrix>0,"+","-"), notecex = 0.75, notecol = "black")
dev.off()

ordered_mrs_disabled <- rev(rownames(final_common_mrs_disabled_activity_matrix)[p2$rowInd])
ordered_p_values_inv_disabled <- NULL
for (i in 1:length(ordered_mrs_disabled))
{
  mr <- ordered_mrs_disabled[i]
  adj_pval <- inv_disabled_median_comparison[inv_disabled_median_comparison$MR==mr,]$Pval
  temp <- cbind(mr,adj_pval)
  ordered_p_values_inv_disabled <- rbind(ordered_p_values_inv_disabled,temp)
}
ordered_p_values_inv_disabled <- as.data.frame(ordered_p_values_inv_disabled)
colnames(ordered_p_values_inv_disabled) <- c("MR","Adj_Pval")
ordered_p_values_inv_disabled$MR <- as.character(as.vector(ordered_p_values_inv_disabled$MR))
ordered_p_values_inv_disabled$Adj.Pval <- as.numeric(as.vector(ordered_p_values_inv_disabled$Adj_Pval))
write.table(ordered_p_values_inv_disabled,"../Results/Paper_Text/Ordered_Pvalues_INV_Disabled_MRs.csv",row.names=F,col.names=T,quote=F,sep=",")

#Get activity of mrs of interest from INV Disabled cancers
mrs_of_interest_disabled <- inv_disabled_median_comparison[inv_disabled_median_comparison$FC_Median<0 &
                                                             inv_disabled_median_comparison$Pval<=0.05,]$MR
disabled_activity_df <- NULL
for (cancer in INV_disabled)
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
  

  for (topmr in mrs_of_interest_disabled)
  {
    high_activity <- amat[topmr,high_indices]
    low_activity <- amat[topmr,low_indices]
    disabled_activity_df <- rbind(disabled_activity_df,
                                  cbind(rep(cancer,length(high_activity)+length(low_activity)),
                                        rep(topmr,length(high_activity)+length(low_activity)),
                                        c(rep("INV High",length(high_activity)),rep("INV Low",length(low_activity))),
                                        c(high_activity,low_activity)))
  }
}

disabled_activity_df <- as.data.frame(disabled_activity_df)
colnames(disabled_activity_df) <- c("Cancer","TopMR","Phenotype","Activity_Value")
disabled_activity_df$Cancer <- as.character(as.vector(disabled_activity_df$Cancer))
disabled_activity_df$TopMR <- as.character(as.vector(disabled_activity_df$TopMR))
disabled_activity_df$Phenotype <- as.character(as.vector(disabled_activity_df$Phenotype))
disabled_activity_df$Activity_Value <- as.numeric(as.vector(disabled_activity_df$Activity_Value))

supp_p4 <- ggplot(data = disabled_activity_df, aes(x=Cancer, y=Activity_Value)) + 
  geom_boxplot(aes(fill=Phenotype)) + facet_wrap( ~ TopMR, nrow=13, ncol=12) + xlab("Cancer") + ylab("Activity Value") +
  geom_point(aes(y=Activity_Value, group=Phenotype), size=0.01, position = position_dodge(width=0.75))+
  guides(fill=guide_legend(title="Phenotype")) + scale_fill_manual(values=c("yellow","green")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.text.x = element_text(angle = -90)) +
  ggtitle("Activities of Top MRs of specific to INV Low in INV Enabled Cancers (>=2 out of 4)") + theme(text = element_text(size=12)) + theme(plot.title = element_text(hjust = 0.5))
ggsave(file="../Results/Paper_Figures/Svgs_Jpgs/Supp_Figure_5.jpg",plot = supp_p4, device = jpeg(), width = 14, height=12, units = "in", dpi = 300)
dev.off()

inv_disabled_median_comparison <- inv_disabled_median_comparison[order(inv_disabled_median_comparison$FC_Median,decreasing = T),]
final_inv_disabled_median_comparison <- inv_disabled_median_comparison[inv_disabled_median_comparison$Pval<=0.05,]
write.table(final_inv_disabled_median_comparison,file="../Results/Paper_Text/All_MRs_TFs_INV_Disabled_Supplementary_Table_S7.csv",
            row.names = F, col.names = T, sep = " & ", quote=F)
#=========================================================================================================================


#Activity of MRs specific to ICR High in all 12 cancers of interest (Figures 4C, 4E)
###########################################################################################################################
mrs_inv_high_all_cancers <- union(final_inv_enabled_median_comparison[final_inv_enabled_median_comparison$FC_Median>0,]$MR,
                                  final_inv_disabled_median_comparison[final_inv_disabled_median_comparison$FC_Median>0,]$MR)

all_cancers <- c(INV_enabled,INV_disabled)
mrs_inv_high_activity_matrix <- matrix(0,nrow=length(mrs_inv_high_all_cancers),ncol=(length(INV_enabled)+length(INV_disabled)))
rownames(mrs_inv_high_activity_matrix) <- mrs_inv_high_all_cancers
colnames(mrs_inv_high_activity_matrix) <- all_cancers
mrs_high_df <- NULL
for (i in 1:length(mrs_inv_high_all_cancers))
{
  topmr <- mrs_inv_high_all_cancers[i]
  print(paste0("Stuck at MR: ",topmr," at iteration ",i))
  inv_enabled_high_activity_values <- NULL
  inv_disabled_high_activity_values <- NULL
  for (j in 1:length(all_cancers))
  {
    cancer <- all_cancers[j]
    if (cancer %in% INV_enabled)
    {
      out_activity <- get_tf_activity_high_low(cancer,topmr)
      inv_enabled_high_activity_values <- c(inv_enabled_high_activity_values,out_activity[[1]])
      mrs_inv_high_activity_matrix[topmr,cancer] <- median(out_activity[[1]])
    } 
    else if (cancer %in% INV_disabled)
    {
      out_activity <- get_tf_activity_high_low(cancer,topmr)
      inv_disabled_high_activity_values <- c(inv_disabled_high_activity_values,out_activity[[1]])
      mrs_inv_high_activity_matrix[topmr,cancer] <- median(out_activity[[1]])
    }
  }
  median_inv_enabled_high <- median(inv_enabled_high_activity_values)
  median_inv_disabled_high <- median(inv_disabled_high_activity_values)
  fc_info <- median(inv_enabled_high_activity_values)-median(inv_disabled_high_activity_values)
  pval <- wilcox.test(inv_enabled_high_activity_values,inv_disabled_high_activity_values,exact=F)$p.value
  temp <- cbind(topmr,fc_info,median_inv_enabled_high,median_inv_disabled_high,pval)
  mrs_high_df <- rbind(mrs_high_df,temp)
}
mrs_high_df <- as.data.frame(mrs_high_df)
colnames(mrs_high_df) <- c("MR","FC_Median","Median_INV_Enabled","Median_INV_Disabled","Padj")
mrs_high_df$MR <- as.character(as.vector(mrs_high_df$MR))
mrs_high_df$Padj <- p.adjust(as.numeric(as.vector(mrs_high_df$Padj)),method="fdr")
mrs_high_df$FC_Median <- round(as.numeric(as.vector(mrs_high_df$FC_Median)),3)
mrs_high_df$Median_INV_Enabled <- round(as.numeric(as.vector(mrs_high_df$Median_INV_Enabled)),3)
mrs_high_df$Median_INV_Disabled <- round(as.numeric(as.vector(mrs_high_df$Median_INV_Disabled)),3)

top_positive_common_mrs_high <- mrs_high_df[mrs_high_df$Median_INV_Enabled>=0 & mrs_high_df$Median_INV_Disabled>=0,]$MR
top_negative_enabled_positive_disabled_mrs_high <- mrs_high_df[mrs_high_df$Median_INV_Enabled<0.0 & 
                                                                 mrs_high_df$Median_INV_Disabled>0 &
                                                                 mrs_high_df$Padj<=0.05,]$MR


interesting_mr_high <- c(top_positive_common_mrs_high,top_negative_enabled_positive_disabled_mrs_high)
interesting_mr_high_df <- mrs_high_df[mrs_high_df$MR %in% interesting_mr_high,]
interesting_mr_high_df <- interesting_mr_high_df[order(interesting_mr_high_df$Median_INV_Enabled,decreasing=T),]
write.table(interesting_mr_high_df,"../Results/Paper_Text/All_MRS_INV_High_Supplementary_Table_S7c.csv",row.names=F,col.names=T,sep="&", quote=F)


positive_mrs_inv_high_activity_matrix <- mrs_inv_high_activity_matrix[top_positive_common_mrs_high,]
pdne_mrs_inv_high_activity_matrix <- mrs_inv_high_activity_matrix[top_negative_enabled_positive_disabled_mrs_high,] 
imp_mr_high_for_overexpression_analysis <- rownames(positive_mrs_inv_high_activity_matrix)
write.table(imp_mr_high_for_overexpression_analysis,"../Results/Paper_Text/All_MRS_INV_High_Important.csv",row.names=F,col.names=F,quote=F)


colcol <- matrix(0,nrow=length(all_cancers),ncol=1)
colcol[c(1:length(INV_enabled)),1] <- "#FDB100"
colcol[c((length(INV_enabled)+1):length(all_cancers)),1] <- "#660066"

#Figure 4C
#====================================================================================
pdf("../Results/Paper_Figures/Figure4/Common_TopMRs_Activity_INV_High_Figure_4C.pdf",height = 10, width=14, pointsize = 12)
p3 <- heatmap.3(positive_mrs_inv_high_activity_matrix, Rowv = FALSE, Colv=FALSE, col = bluered(100), scale="none", main= "Median Activity of MRs specific to INV High Phenotype for 14 cancers in INV High samples",
                dendrogram = "none", key = TRUE, density.info = "none", KeyValueName = "Activity Value", ColSideColors = colcol, ColSideColorsSize = 2,
                margins = c(6,6), useRaster = FALSE, cexRow = 0.5, cexCol = 1.5, cellnote = ifelse(positive_mrs_inv_high_activity_matrix>=0,"+","-"), notecex = 1.0, notecol = "black")
dev.off()

hc.rows <- hclust(dist(pdne_mrs_inv_high_activity_matrix),method='ward.D2')
colcol <- matrix(0,nrow=length(all_cancers),ncol=2)
colnames(colcol) <- c("Enabled/Disabled","INV High")
colcol[c(1:length(INV_enabled)),1] <- "#FDB100"
colcol[c((length(INV_enabled)+1):length(all_cancers)),1] <- "#660066"
colcol[c(1:14),2] <- "yellow"

#Figure 4E
#====================================================================================
pdf("../Results/Paper_Figures/Figure5/Different_TopMRs_Activity_INV_High_Figure_5A.pdf",height = 10, width=12, pointsize = 12)
par(bg="white")
par(fg="black",col.axis="black",col.main="black",col.lab="black", cex.main=1.0)
heatmap.3(pdne_mrs_inv_high_activity_matrix[hc.rows$order,], Rowv = FALSE, Colv=FALSE, col = bluered(100), scale="none", main= "Median Activity of INV High MRs different between INV Enabled & INV Disabled cancers",
          dendrogram = "none", key = TRUE, density.info = "none", KeyValueName = "Activity Value", ColSideColors = colcol, ColSideColorsSize = 3,
          margins = c(6,6), useRaster = FALSE, cexRow = 1.5, cexCol = 1.5, cellnote = ifelse(pdne_mrs_inv_high_activity_matrix[hc.rows$order,]<=0.0,"-","+"), notecex = 2, notecol = "black")
dev.off()
#==================================================================================================================================

#Activity of MRs specific to INV Low in all 14 cancers of interest (Figures 4D, 4F)
#==============================================================================================================================
mrs_inv_low_all_cancers <- union(final_inv_enabled_median_comparison[final_inv_enabled_median_comparison$FC_Median<0,]$MR,
                                 final_inv_disabled_median_comparison[final_inv_disabled_median_comparison$FC_Median<0,]$MR)

all_cancers <- c(INV_enabled,INV_disabled)
mrs_inv_low_activity_matrix <- matrix(0,nrow=length(mrs_inv_low_all_cancers),ncol=(length(INV_enabled)+length(INV_disabled)))
rownames(mrs_inv_low_activity_matrix) <- mrs_inv_low_all_cancers
colnames(mrs_inv_low_activity_matrix) <- all_cancers
mrs_low_df <- NULL
for (i in 1:length(mrs_inv_low_all_cancers))
{
  topmr <- mrs_inv_low_all_cancers[i]
  inv_enabled_low_activity_values <- NULL
  inv_disabled_low_activity_values <- NULL
  for (j in 1:length(all_cancers))
  {
    cancer <- all_cancers[j]
    if (cancer %in% INV_enabled)
    {
      out_activity <- get_tf_activity_high_low(cancer,topmr)
      inv_enabled_low_activity_values <- c(inv_enabled_low_activity_values,out_activity[[2]])
      mrs_inv_low_activity_matrix[topmr,cancer] <- median(out_activity[[2]])
    } 
    else if (cancer %in% INV_disabled)
    {
      out_activity <- get_tf_activity_high_low(cancer,topmr)
      inv_disabled_low_activity_values <- c(inv_disabled_low_activity_values,out_activity[[2]])
      mrs_inv_low_activity_matrix[topmr,cancer] <- median(out_activity[[2]])
    }
  }
  median_inv_enabled_low <- median(inv_enabled_low_activity_values)
  median_inv_disabled_low <- median(inv_disabled_low_activity_values)
  fc_info <- median(inv_enabled_low_activity_values)-median(inv_disabled_low_activity_values)
  pval <- wilcox.test(inv_enabled_low_activity_values,inv_disabled_low_activity_values,exact=F)$p.value
  temp <- cbind(topmr,fc_info,median_inv_enabled_low,median_inv_disabled_low,pval)
  mrs_low_df <- rbind(mrs_low_df,temp)
}
mrs_low_df <- as.data.frame(mrs_low_df)
colnames(mrs_low_df) <- c("MR","FC_Median","Median_INV_Enabled","Median_INV_Disabled","Padj")
mrs_low_df$MR <- as.character(as.vector(mrs_low_df$MR))
mrs_low_df$Padj <- p.adjust(as.numeric(as.vector(mrs_low_df$Padj)),method="fdr")
mrs_low_df$FC_Median <- round(as.numeric(as.vector(mrs_low_df$FC_Median)),3)
mrs_low_df$Median_INV_Enabled <- round(as.numeric(as.vector(mrs_low_df$Median_INV_Enabled)),3)
mrs_low_df$Median_INV_Disabled <- round(as.numeric(as.vector(mrs_low_df$Median_INV_Disabled)),3)

top_positive_common_mrs_low <- mrs_low_df[mrs_low_df$Median_INV_Enabled>0 & mrs_low_df$Median_INV_Disabled>0,]$MR
top_positive_enabled_negative_disabled_mrs_low <- mrs_low_df[mrs_low_df$Median_INV_Enabled>=0.0 & 
                                                               mrs_low_df$Median_INV_Disabled<0 &
                                                               mrs_low_df$Padj<=0.05,]$MR

interesting_mr_low <- c(top_positive_common_mrs_low,top_positive_enabled_negative_disabled_mrs_low)
interesting_mr_low_df <- mrs_low_df[mrs_low_df$MR %in% interesting_mr_low,]
interesting_mr_low_df <- interesting_mr_low_df[order(interesting_mr_low_df$Median_INV_Disabled,decreasing=T),]
write.table(interesting_mr_low_df,"../Results/Paper_Text/All_MRS_INV_Low_Supplementary_Table_S7d.csv",row.names=F,col.names=T,sep="&", quote=F)

positive_mrs_inv_low_activity_matrix <- mrs_inv_low_activity_matrix[top_positive_common_mrs_low,]
nepd_mrs_inv_low_activity_matrix <- mrs_inv_low_activity_matrix[top_positive_enabled_negative_disabled_mrs_low,]
imp_mr_low_for_overexpression_analysis <- rownames(positive_mrs_inv_low_activity_matrix)
write.table(imp_mr_low_for_overexpression_analysis,"../Results/Paper_Text/All_MRS_INV_Low_Important.csv",row.names=F, col.names=F, quote=F)


colcol <- matrix(0,nrow=length(all_cancers),ncol=1)
colcol[c(1:length(INV_enabled)),1] <- "#FDB100"
colcol[c((length(INV_enabled)+1):length(all_cancers)),1] <- "#660066"

#Figure 4D
#====================================================================================
pdf("../Results/Paper_Figures/Figure4/Common_TopMRs_Activity_INV_Low_Figure_4D.pdf",height = 10, width=14, pointsize = 12)
p3 <- heatmap.3(positive_mrs_inv_low_activity_matrix, Rowv = FALSE, Colv=FALSE, col = bluered(100), scale="none", main= "Median Activity of MRs specific to INV Low Phenotype for 12 cancers in INV High samples",
                dendrogram = "none", key = TRUE, density.info = "none", KeyValueName = "Activity Value", ColSideColors = colcol, ColSideColorsSize = 2,
                margins = c(6,6), useRaster = FALSE, cexRow = 0.65, cexCol = 1.5, cellnote = ifelse(positive_mrs_inv_low_activity_matrix>=0,"+","-"), notecex = 1.0, notecol = "black")
dev.off()

#Plot 5B
#######################################################################################
hc.rows <- hclust(dist(nepd_mrs_inv_low_activity_matrix),method='ward.D2')
colcol <- matrix(0,nrow=length(all_cancers),ncol=2)
colnames(colcol) <- c("Enabled/Disabled","INV High")
colcol[c(1:length(INV_enabled)),1] <- "#FDB100"
colcol[c((length(INV_enabled)+1):length(all_cancers)),1] <- "#660066"
colcol[c(1:14),2] <- "green"

pdf("../Results/Paper_Figures/Figure5/Different_TopMRs_Activity_INV_Low_Figure_5B.pdf",height = 10, width=12, pointsize = 12)
par(bg="white")
par(fg="black",col.axis="black",col.main="black",col.lab="black", cex.main=1.0)
heatmap.3(nepd_mrs_inv_low_activity_matrix[hc.rows$order,], Rowv = FALSE, Colv=FALSE, col = bluered(100), scale="none", main= "Median Activity of INV Low MRs different between INV Enabled & INV Disabled cancers",
          dendrogram = "none", key = TRUE, density.info = "none", KeyValueName = "Activity Value", ColSideColors = colcol, ColSideColorsSize = 3,
          margins = c(6,6), useRaster = FALSE, cexRow = 1.5, cexCol = 1.5, cellnote = ifelse(nepd_mrs_inv_low_activity_matrix[hc.rows$order,]<0.0,"-","+"), notecex = 2, notecol = "black")
dev.off()

#########################################################################################################################


#Make figures for downstream enrichment analysis (Supp 6A, 6B)
############################################################################################################
inv_high_goterms_df <- read.table("../Results/Paper_Text/Final_Enriched_GOTerms_INV_High.csv",header=TRUE,sep="\t")
inv_high_goterms_df$generatio <- unlist(lapply(strsplit(as.character(inv_high_goterms_df$members_input_overlap_geneids),split="; "), length))/inv_high_goterms_df$size
temp_df <- inv_high_goterms_df[-log10(inv_high_goterms_df$p.value)>50,]

#Supplementary Figure S7a
supp_p7a <- ggplot(inv_high_goterms_df,aes(x=100*generatio,y=-log10(p.value), shape=term_category)) + geom_point(aes(color=as.factor(term_level))) + 
  xlab("Percentage of MRs in each GO Term") + ylab("-log10(Pvalue)") + scale_shape_manual(name = "GO Terms", values = 1:nlevels(as.factor(inv_high_goterms_df$term_category))) +
  scale_color_manual(name = "Term Level", values =  c("grey","violet","orange")) +
  geom_hline(yintercept=c(0,25,50), color=c("black","blue","red")) + geom_vline(xintercept = 0) + 
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlim(c(0.0,100))+
  geom_text(aes(x=100*generatio,y=-log10(p.value),label=term_name),
            data=temp_df,
            hjust = 0, vjust = 0, nudge_x = 0.025, angle = 0,
            size = 3, check_overlap = T) +
  ggtitle("Enriched GO Terms for INV High Phenotype") + theme(text = element_text(size=14)) + theme(plot.title = element_text(hjust = 0.5))
ggsave(file="../Results/Paper_Figures/Svgs_Jpgs/Supp_Figure_6A.jpg",plot = supp_p7a, device = jpeg(), width = 10, height=10, units = "in", dpi = 300)
dev.off()

##################
inv_low_goterms_df <- read.table("../Results/Paper_Text/Final_Enriched_GOTerms_INV_Low.csv",header=TRUE,sep="\t")
inv_low_goterms_df$generatio <- unlist(lapply(strsplit(as.character(inv_low_goterms_df$members_input_overlap_geneids),split="; "), length))/inv_low_goterms_df$size
temp_df <- inv_low_goterms_df[-log10(inv_low_goterms_df$p.value)>30,]

#Supplementary Figure S7b
supp_p7b <- ggplot(inv_low_goterms_df,aes(x=100*generatio,y=-log10(p.value), shape=term_category)) + geom_point(aes(color=as.factor(term_level))) + 
  xlab("Percentage of MRs in each GO Term") + ylab("-log10(Pvalue)") + scale_shape_manual(name = "GO Terms", values = 1:nlevels(as.factor(inv_low_goterms_df$term_category))) +
  scale_color_manual(name = "Term Level", values =  c("grey","violet","orange")) +
  geom_hline(yintercept=c(0,25,50), color=c("black","blue","red")) + geom_vline(xintercept = 0) + 
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text(aes(x=100*generatio,y=-log10(p.value),label=term_name),
            data=temp_df,
            hjust = 0, vjust = 0, nudge_x = 0.025, angle = 0,
            size = 3, check_overlap = T) +
  coord_cartesian(xlim =c(0, 80), ylim = c(0, 70)) +
  ggtitle("Enriched GO Terms for INV Low Phenotype") + theme(text = element_text(size=14)) + theme(plot.title = element_text(hjust = 0.5))
ggsave(file="../Results/Paper_Figures/Svgs_Jpgs/Supp_Figure_6B.jpg",plot = supp_p7b, device = jpeg(), width = 10, height=10, units = "in", dpi = 300)
dev.off()

#Make Figure 6A_2
#===========================================================================================================
colors <- c('#3cb44b', '#4363d8', '#f58231', '#911eb4', '#f032e6', '#46f0f0' , '#bcf60c', '#a9a9a9', '#e6194b', '#ffe119', '#ffe119', '#e6beff')
inv_low_pathways_df <- read.table("../Results/Paper_Text/Final_Enriched_Pathways_INV_Low.csv",header=TRUE,sep="\t")
inv_low_pathways_df$Description <- as.character(as.vector(inv_low_pathways_df$Description))
inv_low_pathways_df$GeneRatio <- as.character(as.vector(inv_low_pathways_df$GeneRatio))
inv_low_pathways_df$genes <- as.character(as.vector(inv_low_pathways_df$genes))
inv_low_pathways_df$Description <- paste0(inv_low_pathways_df$Description," [",inv_low_pathways_df$GeneRatio,"] ")
inv_low_pathways_df <- inv_low_pathways_df[order(inv_low_pathways_df$generatio),]
p5 <- ggplot(data=inv_low_pathways_df,aes(x=generatio,y=reorder(Description,generatio),size=-log10(pvalue)))+
  geom_point(aes(color=-log10(pvalue)))+xlab("Gene Ratio") + ylab("Enriched Pathways") +  scale_color_continuous(name="-log10(P.adjust)", low="blue", high="red", guide=guide_colorbar(reverse=TRUE))+
  scale_size_continuous(name = "-log10(P.adjust)")+   guides(color=guide_legend(), size = guide_legend())+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background=element_rect(fill = "black"),axis.line = element_line(colour = "white"),
  #  panel.background = element_rect(fill = 'black'),legend.background = element_rect(fill = "black", color = NA),
  #  legend.key = element_rect(color = "gray", fill = "black"), legend.title = element_text(color = "white"), legend.text = element_text(color = "white"))+
  ggtitle("Enriched Pathways for INV Low Phenotype") + theme(text = element_text(size=16,color="black")) + theme(axis.text = element_text(size=14, color="black")) + theme(plot.title = element_text(hjust = 0.5,color="black")) +
  theme(axis.text.y=element_text(color=colors[(inv_low_pathways_df$clusters)]))
ggsave(filename = "../Results/Paper_Figures/Figure6/Enriched_Pathways_INV_Low_Figure_6A_2.pdf",plot=p5,device=pdf(), height=12, width=12, units="in", dpi=300)
dev.off()

#Supplementary Figure S7A
#================================================================================================================
inv_high_pathways_df <- read.table("../Results/Paper_Text/Final_Enriched_Pathways_INV_High.csv",header=TRUE,sep="\t")
inv_high_pathways_df$Description <- as.character(as.vector(inv_high_pathways_df$Description))
inv_high_pathways_df$GeneRatio <- as.character(as.vector(inv_high_pathways_df$GeneRatio))
inv_high_pathways_df$Description <- paste0(inv_high_pathways_df$Description," [",inv_high_pathways_df$GeneRatio,"] ")
inv_high_pathways_df$genes <- as.character(as.vector(inv_high_pathways_df$genes))
inv_high_pathways_df <- inv_high_pathways_df[order(inv_high_pathways_df$generatio),]
supp_p8a <- ggplot(data=inv_high_pathways_df,aes(x=generatio,y=reorder(Description,generatio),size=-log10(pvalue)))+
  geom_point(aes(color=-log10(pvalue)))+xlab("Gene Ratio") + ylab("Enriched Pathways") +  scale_color_continuous(name="-log10(P.adjust)", low="blue", high="red", guide=guide_colorbar(reverse=TRUE))+
  scale_size_continuous(name = "-log10(P.adjust)")+ guides(color=guide_legend(), size = guide_legend())+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Enriched Pathways for INV High Phenotype") + theme(text = element_text(size=10)) + theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y=element_text(color=colors[as.factor(inv_high_pathways_df$clusters)]))
ggsave(filename = "../Results/Paper_Figures/Svgs_Jpgs/Supp_Figure_7A.jpg",plot=supp_p8a,device=jpeg(), height=20, width=12, units="in", dpi=300)
dev.off()

#Supplementary Figure S8B
#=====================================================================================================================
library(circlize)
col_fun = colorRamp2(c(-0.1, 0, 0.1), c("blue", "white", "red"))
col_fun(seq(-3, 3))

inv_high_pathways_df <- inv_high_pathways_df[order(inv_high_pathways_df$generatio,decreasing = T),]
inv_high_pathways_gene_involved <- matrix(0,nrow=nrow(inv_high_pathways_df),ncol=length(rownames(positive_mrs_inv_high_activity_matrix)))
rownames(inv_high_pathways_gene_involved) <- inv_high_pathways_df$Description
colnames(inv_high_pathways_gene_involved) <- rownames(positive_mrs_inv_high_activity_matrix)

for (i in 1:nrow(inv_high_pathways_df))
{
  genes_involved <- unlist(strsplit(inv_high_pathways_df[i,]$genes,split="; "))
  median_activity_scores <- rowMedians(positive_mrs_inv_high_activity_matrix[genes_involved,])
  inv_high_pathways_gene_involved[i,genes_involved] <- median_activity_scores
}
inv_high_pathways_gene_involved <- inv_high_pathways_gene_involved[,which(colSums(inv_high_pathways_gene_involved)>0)]

jpeg("../Results/Paper_Figures/Svgs_Jpgs/Supp_Figure_7B.jpg",height = 900, width=1500, units="px",pointsize = 12)
Heatmap(inv_high_pathways_gene_involved,cluster_rows = FALSE, cluster_columns = FALSE, col=col_fun, row_names_side = "left", 
        name="Median Activity", row_names_gp = gpar(fontsize = 8, col=colors[as.factor(inv_high_pathways_df$clusters)]), column_names_gp = gpar(fontsize = 8), 
        column_title = "Median Activity of MRs in Enriched Pathways specific to INV High")
dev.off()

#Figure 6A_1
#================================================================================================================
rev_inv_low_pathways_df <- inv_low_pathways_df[order(inv_low_pathways_df$generatio,decreasing = T),]
inv_low_pathways_gene_involved <- matrix(0,nrow=nrow(rev_inv_low_pathways_df),ncol=length(rownames(positive_mrs_inv_low_activity_matrix)))
rownames(inv_low_pathways_gene_involved) <- rev_inv_low_pathways_df$Description
colnames(inv_low_pathways_gene_involved) <- rownames(positive_mrs_inv_low_activity_matrix)
for (i in 1:nrow(rev_inv_low_pathways_df))
{
  genes_involved <- unlist(strsplit(rev_inv_low_pathways_df[i,]$genes,split="; "))
  median_activity_scores <- rowMedians(positive_mrs_inv_low_activity_matrix[genes_involved,])
  inv_low_pathways_gene_involved[i,genes_involved] <- median_activity_scores
}

#Make the Sankey plot
inv_low_pathways_gene_involved_to_consider <- inv_low_pathways_gene_involved[,which(colSums(inv_low_pathways_gene_involved)>0)]

#Figure S7C
################################################################################################
jpeg("../Results/Paper_Figures/Svgs_Jpgs/Supp_Figure_7C.jpg",height = 800, width=1200, units="px",pointsize = 12, res = 0.75)
op <- par(oma=c(10,7,1,1))
Heatmap(inv_low_pathways_gene_involved_to_consider,cluster_rows = FALSE, cluster_columns = FALSE, col=col_fun, row_names_side = "left", 
        name="Activity", row_names_gp = gpar(fontsize = 8, col=colors[as.factor(rev_inv_low_pathways_df$clusters)]), column_names_gp = gpar(fontsize = 10), 
        column_title = "Median Activity of MRs in Enriched Pathways specific to INV Low")
dev.off()

INV_Low <- list()
INV_Low$nodes <- data.frame(name = c(colnames(inv_low_pathways_gene_involved_to_consider),rownames(inv_low_pathways_gene_involved_to_consider)))
nodesgroup <- c(rep("white",length(colnames(inv_low_pathways_gene_involved_to_consider))),colors[rev_inv_low_pathways_df$clusters])
INV_Low$nodes$nodesgroup <- nodesgroup
edgelist <- NULL
for (i in 1:nrow(inv_low_pathways_gene_involved_to_consider))
{
  for (j in 1:ncol(inv_low_pathways_gene_involved_to_consider))
  {
    pathway <- rownames(inv_low_pathways_gene_involved_to_consider)[i]
    mr <- colnames(inv_low_pathways_gene_involved_to_consider)[j]
    if (inv_low_pathways_gene_involved_to_consider[pathway,mr]>0)
    {
      mr_id <- which(INV_Low$nodes$name==mr)-1
      pathway_id <- which(INV_Low$nodes$name==pathway)-1
      temp <- cbind(mr_id,pathway_id,abs(inv_low_pathways_gene_involved_to_consider[pathway,mr]),INV_Low$nodes[INV_Low$nodes$name==pathway,"nodesgroup"])
      edgelist <- rbind(edgelist,temp)
    }
  }
}
edgelist <- as.data.frame(edgelist)
colnames(edgelist) <- c("source","target","value","type")
edgelist$source <- as.numeric(as.vector(edgelist$source))
edgelist$target <- as.numeric(as.vector(edgelist$target))
edgelist$value <- as.numeric(as.vector(edgelist$value))
edgelist$type <- as.factor(as.vector(edgelist$type))
INV_Low$links <- edgelist

# putting in a data.frame might help see problems
color_scale <- data.frame(
  range = c(rep("white",length(colnames(inv_low_pathways_gene_involved_to_consider))),colors[rev_inv_low_pathways_df$clusters]),
  domain = INV_Low$nodes$name,
  nodes = INV_Low$nodes,
  stringsAsFactors = FALSE
)

p_sankey <- sankeyNetwork(Links = INV_Low$links, Nodes = INV_Low$nodes, Source = "source",
                   Target = "target", Value = "value", LinkGroup = "type", NodeID = "name", NodeGroup = "nodesgroup",
                   units = "", fontSize = 14, nodeWidth = 25, iterations=0, fontFamily = "Arial", 
                   colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);")) 
#widget2png(p, "../Results/Paper_Figures/Svgs_Jpgs/Sankey_Plot_INV_Low_Figure_6A_1.png")



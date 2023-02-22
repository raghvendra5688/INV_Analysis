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
library(func2vis)
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

INV_enabled <- c("LGG","KIRP","PAAD","MESO","KIRC",
                 "COAD","BLCA","STAD","LUAD","OV")
tfs_enabled_activity_info <- NULL

#Get the list of TFs which are common to all the 10 INV related cancers and their median activity in INV High vs INV Low samples per cancer
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

#Select those transcription regulators which are TFs with regulons of size >= 10 in all the INV Enabled Cancers
common_tfs_enabled_list <- names(which(table(tfs_enabled_activity_info$TFs)>=10))
common_tfs_enabled_activity_info <- tfs_enabled_activity_info[tfs_enabled_activity_info$TFs %in% common_tfs_enabled_list,]

#Make a list of how many times does a common MR appear for the INV Enabled cancers
first_cancer <- INV_enabled[1];
list_common_mrs_enabled <- get_common_mrs(first_cancer)
topMRs_enabled <- rep(list("LGG"),length(list_common_mrs_enabled))
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
all_inv_cancers <- unlist(lapply(topMRs_enabled,function(x) paste(x,collapse=" ")))
no_inv_cancers <- unlist(lapply(topMRs_enabled,function(x) length(x)))
enabled_df <- cbind(names(all_inv_cancers),as.character(all_inv_cancers),as.numeric(no_inv_cancers))
enabled_df <- as.data.frame(enabled_df)
colnames(enabled_df) <- c("MR","List_INV_Cancers","No_INV_Cancers")
enabled_df$MR <- as.character(as.vector(enabled_df$MR))
enabled_df$List_INV_Cancers <- as.character(as.vector(enabled_df$List_INV_Cancers))
enabled_df$No_INV_Cancers <- as.character(as.vector(enabled_df$No_INV_Cancers))
enabled_df <- enabled_df[order(enabled_df$No_INV_Cancers,decreasing=T),]

#Selct MRs which are present in atleast 50% of the INV enabled cancers and are also TF for all INV enabled cancers
temp_enabled_df <- enabled_df[enabled_df$No_INV_Cancers>5,]
temp_enabled_df <- temp_enabled_df[temp_enabled_df$MR %in% common_tfs_enabled_list, ]
write.table(temp_enabled_df,"../Results/Paper_Text/All_MRs_TFs_INV_Supplementary_Table_S3.csv",row.names=F,col.names=T,quote=F,sep=",")

#Get list of common MRs which are all TFs in >5 out of 10 cancers
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
colcol[c(1:length(INV_enabled)),1] <- "red"
colcol[c((length(INV_enabled)+1):(2*length(INV_enabled))),1] <- "blue"

#Perform Wilcox ranksum test or Mann-Whitney test to identify the MRs whose activity between the INV High and INV Low are significant for INV cancers
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
pdf("../Results/Paper_Figures/Figure4/Common_TopMRs_Activity_INV_Figure_4A.pdf",height = 12, width=14, pointsize = 11)
par(bg="white")
par(fg="black",col.axis="black",col.main="black",col.lab="black", cex.main=1.5)
p1 <- heatmap.3(final_common_mrs_enabled_activity_matrix, Rowv = TRUE, Colv=F, col = bluered(100), scale="none", main= "Median Activity of Common MRs in INV Enabled Cancers", # (>=4 out of 8)",
                dendrogram = "both", key = TRUE, density.info = "none", KeyValueName = "Activity Value", ColSideColors = colcol, ColSideColorsSize = 2,
                margins = c(6,6), useRaster = FALSE, cexRow = 1.0, cexCol = 1.5, cellnote = ifelse(final_common_mrs_enabled_activity_matrix>0,"+","-"), notecex = 1, notecol = "black")
dev.off()

#Identify the MRs which are specific to INV Low
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

supp_p5 <- ggplot(data = enabled_activity_df, aes(x=Cancer, y=Activity_Value)) +  
  geom_boxplot(aes(fill=Phenotype)) + facet_wrap( ~ TopMR, nrow=7, ncol=10) + xlab("Cancer") + ylab("Activity Value") +
  geom_point(aes(y=Activity_Value, group=Phenotype), size=0.01, position = position_dodge(width=0.75))+
  guides(fill=guide_legend(title="Phenotype")) + scale_fill_manual(values=c("red","blue")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.text.x = element_text(angle = -90)) +
  ggtitle("Activities of Top MRs of specific to INV Low in INV Enabled Cancers (>=5 out of 10)") + theme(text = element_text(size=12)) + theme(plot.title = element_text(hjust = 0.5))
ggsave(file="../Results/Paper_Figures/Svgs_Jpgs/Supp_Figure_5.jpg",plot = supp_p5, device = jpeg(), width = 14, height=12, units = "in", dpi = 300)
dev.off()

inv_enabled_median_comparison <- inv_enabled_median_comparison[order(inv_enabled_median_comparison$FC_Median,decreasing = T),]
final_inv_enabled_median_comparison <- inv_enabled_median_comparison[inv_enabled_median_comparison$Pval<=0.05,]
final_inv_enabled_median_comparison$Pval <- as.numeric(as.vector(signif(final_inv_enabled_median_comparison$Pval,digits=6)))
write.table(final_inv_enabled_median_comparison,file="../Results/Paper_Text/All_MRs_TFs_INV_Supplementary_Table_S4.csv",
            row.names = F, col.names = T, sep = " & ", quote=F)

inv_low_specific_mrs <- inv_enabled_median_comparison[inv_enabled_median_comparison$FC_Median<0 &
                                                        inv_enabled_median_comparison$Pval<=0.05,]$MR
inv_high_specific_mrs <- inv_enabled_median_comparison[inv_enabled_median_comparison$FC_Median>0 & 
                                                         inv_enabled_median_comparison$Pval<=0.05,]$MR
write.table(inv_low_specific_mrs,"../Results/Paper_Text/All_INV_Low_Specific_MRs.csv", row.names=F, col.names=F, quote=F)
write.table(inv_high_specific_mrs,"../Results/Paper_Text/All_INV_High_Specific_MRs.csv", row.names=F, col.names=F, quote=F)
write.table(rownames(D),"../Results/Paper_Text/Background_genes.csv", row.names=F, col.names=F, quote=F)
########################################################################################################################################################

#Perform downstream GO Term enrichment and pathway enrichment for INV High and INV Low specific MRs

#Make figures for downstream enrichment analysis (Supp 6A, 6B)
############################################################################################################
inv_high_goterms_df <- read.table("../Results/Paper_Text/Enriched_GOTerms_INV_High.csv",header=TRUE,sep="\t")
inv_high_goterms_df$generatio <- unlist(lapply(strsplit(as.character(inv_high_goterms_df$members_input_overlap_geneids),split="; "), length))/inv_high_goterms_df$size
inv_high_case_vs_control <- data.frame(gene=inv_high_specific_mrs,fc=rep(2,length(inv_high_specific_mrs)))
final_inv_high_goterms_df <- clean_go_terms(df_case_vs_ctrl = inv_high_case_vs_control, df_goterms = inv_high_goterms_df)

temp_df <- final_inv_high_goterms_df[-log10(final_inv_high_goterms_df$p.value)>35,]

#Supplementary Figure S6a
supp_p6a <- ggplot(final_inv_high_goterms_df,aes(x=100*generatio,y=-log10(p.value), shape=term_category)) + geom_point(aes(color=as.factor(term_level))) + 
  xlab("Percentage of MRs in each GO Term") + ylab("-log10(Pvalue)") + scale_shape_manual(name = "GO Terms", values = 1:nlevels(as.factor(inv_high_goterms_df$term_category))) +
  scale_color_manual(name = "Term Level", values =  c("grey","violet","orange")) +
  geom_hline(yintercept=c(1,35,50), color=c("black","blue","red")) + geom_vline(xintercept = 0) + 
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlim(c(0.0,100))+
  geom_text(aes(x=100*generatio,y=-log10(p.value),label=term_name),
            data=temp_df,
            hjust = 0, vjust = 0, nudge_x = 0.025, angle = 0,
            size = 3, check_overlap = T) +
  ggtitle("Enriched GO Terms for INV High Phenotype") + theme(text = element_text(size=14)) + theme(plot.title = element_text(hjust = 0.5))
ggsave(file="../Results/Paper_Figures/Svgs_Jpgs/Supp_Figure_6A.pdf",plot = supp_p6a, device = pdf(), width = 10, height=10, units = "in", dpi = 300)
dev.off()

###############################################################################
inv_low_goterms_df <- read.table("../Results/Paper_Text/Enriched_GOTerms_INV_Low.csv",header=TRUE,sep="\t")
inv_low_goterms_df$generatio <- unlist(lapply(strsplit(as.character(inv_low_goterms_df$members_input_overlap_geneids),split="; "), length))/inv_low_goterms_df$size
inv_low_case_vs_control <- data.frame(gene=inv_low_specific_mrs,fc=rep(-1,length(inv_low_specific_mrs)))
final_inv_low_goterms_df <- clean_go_terms(df_case_vs_ctrl = inv_low_case_vs_control, df_goterms = inv_low_goterms_df)
temp_df <- inv_low_goterms_df[-log10(inv_low_goterms_df$p.value)>20,]

#Supplementary Figure S7b
supp_p6b <- ggplot(inv_low_goterms_df,aes(x=100*generatio,y=-log10(p.value), shape=term_category)) + 
  geom_point(aes(color=as.factor(term_level))) + 
  xlab("Percentage of MRs in each GO Term") + ylab("-log10(Pvalue)") + 
  scale_shape_manual(name = "GO Terms", values = 1:nlevels(as.factor(inv_low_goterms_df$term_category))) +
  scale_color_manual(name = "Term Level", values =  c("grey","violet","orange")) +
  geom_hline(yintercept=c(0,20,30), color=c("black","blue","red")) + geom_vline(xintercept = 0) + 
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text(aes(x=100*generatio,y=-log10(p.value),label=term_name),
            data=temp_df,
            hjust = 0, vjust = 0, nudge_x = 0.025, angle = 0,
            size = 3, check_overlap = T) +
  coord_cartesian(xlim =c(0, 80), ylim = c(0, 30)) +
  ggtitle("Enriched GO Terms for INV Low Phenotype") + theme(text = element_text(size=14)) + theme(plot.title = element_text(hjust = 0.5))
ggsave(file="../Results/Paper_Figures/Svgs_Jpgs/Supp_Figure_6B.pdf",plot = supp_p6b, device = pdf(), width = 10, height=10, units = "in", dpi = 300)
dev.off()

#Make Figure 4B
#===========================================================================================================
colors <- c('#3cb44b', '#4363d8', '#f58231', '#911eb4', '#f032e6', '#46f0f0' , '#bcf60c', '#a9a9a9', '#e6194b', '#ffe119', '#ffe119', '#e6beff')
inv_low_pathways_df <- read.table("../Results/Paper_Text/Enriched_Pathways_INV_Low.csv", header=T, sep="\t")
final_inv_low_pathways_df <- clean_pathways(df_case_vs_ctrl = inv_low_case_vs_control, df_pathway = inv_low_pathways_df)
write.table(final_inv_low_pathways_df,"../Results/Paper_Text/Supplementary_Table_S5_Final_Enriched_INV_Low_Pathways.csv",row.names=F, col.names=T, quote=F, sep="\t")

final_inv_low_pathways_df$pathway <- as.character(as.vector(final_inv_low_pathways_df$pathway))
final_inv_low_pathways_df$Generatio <- paste0(final_inv_low_pathways_df$genes_down,"/",final_inv_low_pathways_df$effective_size)
final_inv_low_pathways_df$generatio <- final_inv_low_pathways_df$genes_down/final_inv_low_pathways_df$effective_size
final_inv_low_pathways_df$Description <- paste0(final_inv_low_pathways_df$pathway,' [',final_inv_low_pathways_df$Generatio,']')

final_inv_low_pathways_df <- final_inv_low_pathways_df[order(final_inv_low_pathways_df$generatio),]
p5 <- ggplot(data=final_inv_low_pathways_df,aes(x=generatio,y=reorder(Description,generatio),size=-log10(q.value)))+
  geom_point(aes(color=-log10(q.value)))+xlab("Gene Ratio") + ylab("Enriched Pathways") +  
  scale_color_continuous(name="-log10(P.adjust)", low="blue", high="red", guide=guide_colorbar(reverse=TRUE))+
  scale_size_continuous(name = "-log10(P.adjust)")+   guides(color=guide_legend(), size = guide_legend())+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Enriched Pathways for INV Low Phenotype") + 
  theme(text = element_text(size=14,color="black")) + theme(axis.text = element_text(size=14, color="black")) + 
  theme(plot.title = element_text(hjust = 0.5,color="black")) +
  theme(axis.text.y=element_text(color=colors[(final_inv_low_pathways_df$clusters)]))
ggsave(filename = "../Results/Paper_Figures/Figure4/Enriched_Pathways_INV_Low_Figure_4B.pdf",plot=p5,device=pdf(), height=6, width=10, units="in", dpi=300)
dev.off()

#Supplementary Figure 4C
#================================================================================================================
inv_high_pathways_df <- read.table("../Results/Paper_Text/Enriched_Pathways_INV_High.csv",header=TRUE,sep="\t")
final_inv_high_pathways_df <- clean_pathways(df_case_vs_ctrl = inv_high_case_vs_control, df_pathway = inv_high_pathways_df)
write.table(final_inv_high_pathways_df,"../Results/Paper_Text/Supplementary_Table_S6_Final_Enriched_Pathways_INV_High.csv",row.names=F, col.names=T, quote=F,sep="\t")

final_inv_high_pathways_df$pathway <- as.character(as.vector(final_inv_high_pathways_df$pathway))
final_inv_high_pathways_df$Generatio <- paste0(final_inv_high_pathways_df$genes_up,"/",final_inv_high_pathways_df$effective_size)
final_inv_high_pathways_df$generatio <- final_inv_high_pathways_df$genes_up/final_inv_high_pathways_df$effective_size
final_inv_high_pathways_df$Description <- paste0(final_inv_high_pathways_df$pathway,' [',final_inv_high_pathways_df$Generatio,']')

final_inv_high_pathways_df <- final_inv_high_pathways_df[order(final_inv_high_pathways_df$p.value),]
top_final_inv_high_pathways_df <- final_inv_high_pathways_df[c(1:30),]

top_final_inv_high_pathways_df <- top_final_inv_high_pathways_df[order(top_final_inv_high_pathways_df$generatio),]
p6 <- ggplot(data=top_final_inv_high_pathways_df,aes(x=generatio,y=reorder(Description,generatio),size=-log10(q.value)))+
  geom_point(aes(color=-log10(q.value)))+xlab("Gene Ratio") + ylab("Enriched Pathways") +  
  scale_color_continuous(name="-log10(P.adjust)", low="blue", high="red", guide=guide_colorbar(reverse=TRUE))+
  scale_size_continuous(name = "-log10(P.adjust)")+ guides(color=guide_legend(), size = guide_legend())+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Enriched Pathways for INV High Phenotype") + theme(text = element_text(size=14)) + theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y=element_text(color=colors[as.factor(top_final_inv_high_pathways_df$clusters)]))
ggsave(filename = "../Results/Paper_Figures/Figure4/Enriched_Pathways_INV_High_Figure_4D.pdf",plot=p6,device=jpeg(), height=12, width=12, units="in", dpi=300)
dev.off()

final_inv_high_pathways_df <- final_inv_high_pathways_df[order(final_inv_high_pathways_df$generatio),]
p7 <- ggplot(data=final_inv_high_pathways_df,aes(x=generatio,y=reorder(Description,generatio),size=-log10(q.value)))+
  geom_point(aes(color=-log10(q.value)))+xlab("Gene Ratio") + ylab("Enriched Pathways") +  
  scale_color_continuous(name="-log10(P.adjust)", low="blue", high="red", guide=guide_colorbar(reverse=TRUE))+
  scale_size_continuous(name = "-log10(P.adjust)")+ guides(color=guide_legend(), size = guide_legend())+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Enriched Pathways for INV High Phenotype") + theme(text = element_text(size=14)) + theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y=element_text(color=colors[as.factor(final_inv_high_pathways_df$clusters)]))
ggsave(filename = "../Results/Paper_Figures/Svgs_Jpgs/Supp_Figure_8.pdf",plot=p7,device=pdf(), height=12, width=12, units="in", dpi=300)
dev.off()


#Supplementary Figure S8B
#=====================================================================================================================
library(circlize)
col_fun = colorRamp2(c(-0.1, 0, 0.1), c("blue", "white", "red"))
col_fun(seq(-3, 3))

#Figure 4C
#================================================================================================================
positive_mrs_inv_high_activity_matrix <- common_mrs_enabled_activity_matrix[inv_high_specific_mrs,]
rev_inv_high_pathways_df <- top_final_inv_high_pathways_df[order(top_final_inv_high_pathways_df$generatio,decreasing = T),]
inv_high_pathways_gene_involved <- matrix(0,nrow=nrow(rev_inv_high_pathways_df),ncol=length(rownames(positive_mrs_inv_high_activity_matrix)))
rownames(inv_high_pathways_gene_involved) <- rev_inv_high_pathways_df$Description
colnames(inv_high_pathways_gene_involved) <- rownames(positive_mrs_inv_high_activity_matrix)
for (i in 1:nrow(rev_inv_high_pathways_df))
{
  genes_involved <- unlist(strsplit(rev_inv_high_pathways_df[i,]$members_input_overlap,split="; "))
  median_activity_scores <- rowMedians(positive_mrs_inv_high_activity_matrix[genes_involved,])
  inv_high_pathways_gene_involved[i,genes_involved] <- median_activity_scores
}

#Make the Sankey plot
inv_high_pathways_gene_involved_to_consider <- inv_high_pathways_gene_involved[,which(colSums(inv_high_pathways_gene_involved)>0)]

#Figure S7
################################################################################################
pdf("../Results/Paper_Figures/Svgs_Jpgs/Supp_Figure_7.pdf",height = 8, width=12)
op <- par(oma=c(10,7,1,1))
Heatmap(inv_high_pathways_gene_involved_to_consider,cluster_rows = FALSE, cluster_columns = FALSE, col=col_fun, row_names_side = "left", 
        name="Activity", row_names_gp = gpar(fontsize = 10, col=colors[as.factor(rev_inv_high_pathways_df$clusters)]), column_names_gp = gpar(fontsize = 10), 
        column_title = "Median Activity of MRs in Enriched Pathways specific to INV High")
dev.off()

INV_High <- list()
INV_High$nodes <- data.frame(name = c(colnames(inv_high_pathways_gene_involved_to_consider),rownames(inv_high_pathways_gene_involved_to_consider)))
nodesgroup <- c(rep("white",length(colnames(inv_high_pathways_gene_involved_to_consider))),colors[rev_inv_high_pathways_df$clusters])
INV_High$nodes$nodesgroup <- nodesgroup
edgelist <- NULL
for (i in 1:nrow(inv_high_pathways_gene_involved_to_consider))
{
  for (j in 1:ncol(inv_high_pathways_gene_involved_to_consider))
  {
    pathway <- rownames(inv_high_pathways_gene_involved_to_consider)[i]
    mr <- colnames(inv_high_pathways_gene_involved_to_consider)[j]
    if (inv_high_pathways_gene_involved_to_consider[pathway,mr]>0)
    {
      mr_id <- which(INV_High$nodes$name==mr)-1
      pathway_id <- which(INV_High$nodes$name==pathway)-1
      temp <- cbind(mr_id,pathway_id,abs(inv_high_pathways_gene_involved_to_consider[pathway,mr]),INV_High$nodes[INV_High$nodes$name==pathway,"nodesgroup"])
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
INV_High$links <- edgelist

# putting in a data.frame might help see problems
color_scale <- data.frame(
  range = c(rep("white",length(colnames(inv_high_pathways_gene_involved_to_consider))),colors[rev_inv_high_pathways_df$clusters]),
  domain = INV_High$nodes$name,
  nodes = INV_High$nodes,
  stringsAsFactors = FALSE
)

p_sankey <- sankeyNetwork(Links = INV_High$links, Nodes = INV_High$nodes, Source = "source",
                   Target = "target", Value = "value", LinkGroup = "type", NodeID = "name", NodeGroup = "nodesgroup",
                   units = "", fontSize = 14, nodeWidth = 25, iterations=0, fontFamily = "Arial", 
<<<<<<< HEAD
                   colourScale = JS('d3.scaleOrdinal().domain(["Fristående","Nästan_öppet","Halvöppet"])
                                    .range(["#EDF8E9","#BAE4B3","#74C476"])'))
=======
                   colourScale = JS('d3.scaleOrdinal().domain(["Fristående","Nästan_öppet","Halvöppet","Slutet"])
                                    .range(["#EDF8E9","#BAE4B3","#74C476","#238B45"])'))
>>>>>>> 686c8ab687c42f7129469caff3a2dae3da6e9efa
                   #colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);")) 
#widget2png(p, "../Results/Paper_Figures/Svgs_Jpgs/Sankey_Plot_INV_Low_Figure_6A_1.png")



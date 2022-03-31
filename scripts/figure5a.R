library(data.table)
library(ggplot2)
library(gplots)
library(Matrix)
library(erer)
library(doMC)
library(heatmap3)
library(mgcv)
library(devtools)

#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

setwd('../scripts/')

source('gene-reverse-network.R')
source('get_functions.R')

all_cold_mrs_df <- read.table("../Results/Paper_Text/All_INV_Low_Specific_MRs.csv",header=F)
all_cold_mrs <- as.character(as.vector(all_cold_mrs_df$V1))

all_hot_mrs_df <- read.table("../Results/Paper_Text/All_INV_High_Specific_MRs.csv",header=F)
all_hot_mrs <- as.character(as.vector(all_hot_mrs_df$V1))
all_mrs <- c(all_hot_mrs,all_cold_mrs)

RGBM_path <- c("../Data/");

inv_neutral_cancers <- setdiff(setdiff(list.files(RGBM_path),inv_prognostic_cancers),c("Others","Assembler_Panca_Normalized","ARACNE","PanCancer","PRECOG","survival_data"))


get_activity_matrix_info <- function(counter,sample_ids,cancers,activity_matrix,cancer_info_vector,inv_pheno_info_vector,inv_type_info_vector,inv_type)
{
  for (i in 1:length(cancers))
  {
    cancer_type <- cancers[i]
    load(paste0("../Results/",cancer_type,"/Adjacency_Matrix/",cancer_type,"_Full_Activity_matrix_FGSEA.Rdata"))
    
    #Load mechanistic network
    #======================================================================================
    load('../Data/Others/me_net_full.Rdata')
    
    #Get the data
    filename = cancer_type
    out <- loading_data(filename,M)
    D <- as.matrix(log2(t(out[[1]])+1))
    
    #Get high and low indices
    load(paste0("../Data/",filename,"/",filename,"_INV_cluster_assignment_k2-6_v2.Rdata"))
    high_low_output <- get_high_low_indices(table_cluster_assignment,D)
    table_cluster_assignment <- high_low_output[[1]]
    high_indices <- high_low_output[[2]]
    low_indices <- high_low_output[[3]]
    
    amat_to_use <- amat[,c(high_indices,low_indices)]
    amat_to_use[amat_to_use>0] <- amat_to_use[amat_to_use>0]/max(amat_to_use)
    amat_to_use[amat_to_use<0] <- amat_to_use[amat_to_use<0]/abs(min(amat_to_use))
    samples <- dim(amat_to_use)[2]
    high_index_last <- length(high_indices)
    low_index_last <- ncol(amat_to_use)
    
    cancer_info_vector <- c(cancer_info_vector,rep(cancer_type,samples))
    inv_pheno_info_vector <- c(inv_pheno_info_vector,rep("INV High",length(high_indices)),
                               rep("INV Low",length(low_indices)))
    inv_type_info_vector <- c(inv_type_info_vector,rep(inv_type,samples))
    sample_ids <- c(sample_ids,colnames(amat_to_use))
    for (j in 1:length(all_mrs))
    {
      mr <- all_mrs[j]
      if (mr %in% rownames(amat_to_use))
      {
        activity_matrix[mr,counter+c(1:high_index_last)] <- amat_to_use[mr,c(1:high_index_last)]
        activity_matrix[mr,counter+c((high_index_last+1):low_index_last)] <- amat_to_use[mr,c((high_index_last+1):low_index_last)]
      }
    }
    counter <- counter+samples
  }
  
  output_list <- list()
  output_list[[1]] <- counter
  output_list[[2]] <- cancer_info_vector
  output_list[[3]] <- inv_pheno_info_vector
  output_list[[4]] <- inv_type_info_vector
  output_list[[5]] <- activity_matrix
  output_list[[6]] <- sample_ids
  return(output_list)
}

cancers <- inv_neutral_cancers
cancer_samples <- 0
for (i in 1:length(cancers))
{
  cancer_type <- cancers[i]
  print(paste0("At cancer ",cancer_type))
  #high_indices_table <- read.table(paste0("../Results/",cancer_type,"/Adjacency_Matrix/",cancer_type,"_Full_high_indices.csv"),header=TRUE)
  #low_indices_table <- read.table(paste0("../Results/",cancer_type,"/Adjacency_Matrix/",cancer_type,"_Full_low_indices.csv"),header=TRUE)
  #Load mechanistic network
  #======================================================================================
  load('../Data/Others/me_net_full.Rdata')
  
  #Get the data
  filename = cancer_type
  out <- loading_data(filename,M)
  D <- as.matrix(log2(t(out[[1]])+1))
  
  #Get high and low indices
  load(paste0("../Data/",filename,"/",filename,"_INV_cluster_assignment_k2-6_v2.Rdata"))
  high_low_output <- get_high_low_indices(table_cluster_assignment,D)
  table_cluster_assignment <- high_low_output[[1]]
  high_indices <- high_low_output[[2]]
  low_indices <- high_low_output[[3]]
  
  cancer_samples <- cancer_samples+length(high_indices)+length(low_indices)
}

Activity_Matrix <- Matrix(0,nrow=length(all_mrs),ncol=cancer_samples)
rownames(Activity_Matrix) <- all_mrs
counter=0;
cancer_info_vector <- NULL
inv_pheno_info_vector <- NULL
inv_type_info_vector <- NULL
sample_ids <- NULL

# output <- get_activity_matrix_info(counter,sample_ids,icr_disabled_cancers,Activity_Matrix,cancer_info_vector,icr_pheno_info_vector,icr_type_info_vector,"ICR Disabled")
# counter <- output[[1]]
# cancer_info_vector <- output[[2]]
# icr_pheno_info_vector <- output[[3]]
# icr_type_info_vector <- output[[4]]
# Activity_Matrix <- output[[5]]
# sample_ids <- output[[6]]
# 
# output <- get_activity_matrix_info(counter,sample_ids,icr_enabled_cancers,Activity_Matrix,cancer_info_vector,icr_pheno_info_vector,icr_type_info_vector,"ICR Enabled")
# counter <- output[[1]]
# cancer_info_vector <- output[[2]]
# icr_pheno_info_vector <- output[[3]]
# icr_type_info_vector <- output[[4]]
# Activity_Matrix <- output[[5]]
# sample_ids <- output[[6]]

output <- get_activity_matrix_info(counter,sample_ids,inv_neutral_cancers,Activity_Matrix,cancer_info_vector,inv_pheno_info_vector,inv_type_info_vector,"INV Neutral")
counter <- output[[1]]
cancer_info_vector <- output[[2]]
inv_pheno_info_vector <- output[[3]]
inv_type_info_vector <- output[[4]]
Activity_Matrix <- output[[5]]
sample_ids <- output[[6]]
colnames(Activity_Matrix) <- sample_ids
save(Activity_Matrix,all_cold_mrs,all_hot_mrs,cancer_info_vector,inv_type_info_vector,file="../Results/PanCancer/INV_Neutral_Activity_Info.Rdata")


#Make heatmap using heatmap.3 function
n <- length(unique(cancer_info_vector))
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
set.seed(43)
col=sample(color, n)
Cancers <- col[as.factor(cancer_info_vector)]
inv_colors <- c("red","blue")
INV_Clusters <- inv_colors[as.factor(inv_pheno_info_vector)]
inv_type_colors <- c("#660066","#FDB100")
INV_EorD <- inv_type_colors[as.factor(inv_type_info_vector)]

clab=cbind(Cancers,INV_Clusters,INV_EorD)
rlab = t(c(rep("red",length(all_hot_mrs)),rep("blue",length(all_cold_mrs))))
colnames(clab) <- c("Cancers","INV Phenotype","INV Cancer Clusters")
rownames(rlab) <- "MRs"

#Define custom dist and hclust functions for use with heatmaps
Activity_Matrix <- as.matrix(Activity_Matrix)
#2nd column is INV High vs INV Low
hc.cols1 <- order(clab[,2],decreasing = T)
clab1 <- clab[hc.cols1,]
Activity_Matrix <- Activity_Matrix[,hc.cols1]
newhc.cols <- NULL
#3rd column is INV Prognostic vs INV Neutral clusters
for (i in 1:length(unique(clab[,2])))
{
  color_id <- unique(clab[,2])[i]
  newhc.cols <- c(newhc.cols,length(newhc.cols)+order(clab1[clab1[,2]==color_id,3],decreasing = T))
}
Activity_Matrix <- Activity_Matrix[,newhc.cols]
clab2 <- clab1[newhc.cols,]
lower_val <- quantile(Activity_Matrix,probs=c(0.05))
upper_val <- quantile(Activity_Matrix,probs=c(0.95))
Activity_Matrix[Activity_Matrix<lower_val] <- lower_val
Activity_Matrix[Activity_Matrix>upper_val] <- upper_val

.merge_hclust <- function(hclist) {
  #-- Merge
  d <- as.dendrogram(hclist[[1]])
  for (i in 2:length(hclist)) {
    d <- merge(d, as.dendrogram(hclist[[i]]))
  }
  return(as.hclust(d))
}

all_inv_high_indices <- which(clab2[,2]==inv_colors[1])
all_inv_low_indices <- which(clab2[,2]==inv_colors[2])
hc.cols.inv_high <- hclust(dist(t(Activity_Matrix[,all_inv_high_indices])),method='ward.D2')
hc.cols.inv_low <- hclust(dist(t(Activity_Matrix[,all_inv_low_indices])),method='ward.D2')
hc.cols <- c(hc.cols.inv_high$order,(length(all_inv_high_indices)+hc.cols.inv_low$order))
hc_col_list <- list(hc.cols.inv_high,hc.cols.inv_low)
hc_col <- .merge_hclust(hc_col_list)

all_hot_mr_indices <- which(rlab[1,]==inv_colors[1])
all_cold_mr_indices <- which(rlab[1,]==inv_colors[2])
hc.rows.inv_high <- hclust(dist(Activity_Matrix[all_hot_mr_indices,]),method='ward.D2')
hc.rows.inv_low <- hclust(dist(Activity_Matrix[all_cold_mr_indices,]),method='ward.D2')
hc.rows <- c(hc.rows.inv_high$order,(length(all_hot_mr_indices)+hc.rows.inv_low$order))
hc_row_list <- list(hc.rows.inv_high,hc.rows.inv_low)
hc_row <- .merge_hclust(hc_row_list)


#Main Plotting Function for Figure 5A
################################################################################
pdf("../Results/Paper_Figures/Figure5/Heatmap_All_MRs_Neutral_Cancers_Figure_5A.pdf",height = 13, width=15, pointsize = 14)
par(bg="white")
par(fg="black",col.axis="black",col.main="black",col.lab="black", cex.main=2.0)
p1 <- heatmap.3(Activity_Matrix, scale="none", dendrogram="both", margins=c(4,20),
                Rowv=as.dendrogram(hc_row), Colv=as.dendrogram(hc_col), ColSideColors=clab2, labCol = FALSE,
                key=TRUE, keysize=1.5, labRow=all_mrs, cexRow=0.8, col=colorspace::diverge_hsv(100),RowSideColors=rlab,
                ColSideColorsSize=3, RowSideColorsSize=1, KeyValueName="Activity Value")
legend("topright",legend=c(unique(cancer_info_vector),"","",unique(inv_type_info_vector),"","",unique(inv_pheno_info_vector)),
       fill=c(unique(Cancers),"white","white",unique(INV_EorD),"white","white",unique(INV_Clusters)), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
dev.off()

ordered_mrs <- rev(rownames(Activity_Matrix)[p1$rowInd])
high_ids <- clab2[,2]=="red"
low_ids <- clab2[,2]=="blue"
activity_df <- perform_wilcox_test(Activity_Matrix[ordered_mrs,high_ids],Activity_Matrix[ordered_mrs,low_ids])
activity_df <- activity_df[,c(2:ncol(activity_df))]
activity_df$Mean1 <- round(activity_df$Mean1,3)
activity_df$Mean2 <- round(activity_df$Mean2,3)
activity_df$FC_Mean <- round(activity_df$FC_Mean,3)
activity_df$Pval <- signif(activity_df$Pval,3)
activity_df$Padjust <- signif(activity_df$Padjust,3)
activity_df <- activity_df[order(activity_df$Padjust,decreasing=T),]
write.table(activity_df,"../Results/Paper_Text/Diff_MRS_Activity_INV_Neutral_Table_S7.csv",row.names=T,col.names=T,quote=F,sep=",")

library(data.table)
library(ggplot2)
library(gplots)
library(Matrix)
library(erer)
library(doMC)
library(heatmap3)
library(mgcv)
library(devtools)
library(gprofiler2)
library(Biobase)

#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

setwd('../scripts/')

source('gene-reverse-network.R')
source('get_functions.R')

convert_character_vector <- function(df)
{
  for (i in 1:ncol(df))
  {
    df[,i] <- as.character(as.vector(df[,i]))
  }
  return(df)
}

#For survival analysis
###################################################################################################################
all_cold_mrs_df <- read.table("../Results/Paper_Text/All_INV_Low_Specific_MRs.csv",header=F)
all_cold_mrs <- as.character(as.vector(all_cold_mrs_df$V1))

all_hot_mrs_df <- read.table("../Results/Paper_Text/All_INV_High_Specific_MRs.csv",header=F)
all_hot_mrs <- as.character(as.vector(all_hot_mrs_df$V1))
all_mrs <- c(all_hot_mrs,all_cold_mrs)
RGBM_path <- c("../Data/");

inv_prognostic_cancers <- c("LGG","KIRP","PAAD","MESO","KIRC","COAD","BLCA","STAD","LUAD","OV")

cancers <- inv_prognostic_cancers
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
  valid_mrs <- rownames(amat_to_use)[rownames(amat_to_use) %in% all_mrs]
  inv_pheno_info_vector <- c(rep("INV High",length(high_indices)),
                             rep("INV Low",length(low_indices)))
  sample_ids <- colnames(amat_to_use)
  valid_amat_to_use <- amat_to_use[valid_mrs,]
  save(list = c('valid_amat_to_use','sample_ids','inv_pheno_info_vector'), file=paste0("../Results/PanCancer/",cancer_type,"_activity_info_for_survival.Rdata"))
}

#Analyze the PRECOG dataset
##########################################################################################################
predictions_df <- readRDS("../Data/PRECOG/predictors_all_merged3.rds")
expr_list1 <- readRDS("../Data/PRECOG/es_list1.rds")
expr_list2 <- readRDS("../Data/PRECOG/es_list2.rds")
mapping_df <- read.table("../Data/PRECOG/Mapping_GEO.csv",header=TRUE,sep=",")
mapping_df <- convert_character_vector(mapping_df)

#Get only those samples for which cancer type is available
cancer_list <- NULL
rev_predictions_df <- predictions_df[predictions_df$Study %in% mapping_df$Study,]
for (i in 1:nrow(rev_predictions_df))
{
  study_type <- rev_predictions_df$Study[i]
  cancer_type <- mapping_df[mapping_df$Study==study_type,]$Cancer  
  cancer_list <- c(cancer_list,cancer_type)
}
rev_predictions_df$Cancer <- cancer_list

#Get largest cohort for each cancer type
geo_cancer_mat <- as.matrix(table(rev_predictions_df$Study,rev_predictions_df$Cancer))
unique_cancers <- colnames(geo_cancer_mat)
geo_cancer_edgelist <- NULL
for (i in 1:length(unique_cancers))
{
  id <- which.max(geo_cancer_mat[,colnames(geo_cancer_mat)==unique_cancers[i]])
  geo_id <- rownames(geo_cancer_mat)[id]
  sample_size <- geo_cancer_mat[geo_id,unique_cancers[i]]
  temp <- cbind(geo_id,unique_cancers[i],sample_size)
  geo_cancer_edgelist <- rbind(geo_cancer_edgelist,temp)
}
geo_cancer_edgelist <- as.data.frame(geo_cancer_edgelist)
colnames(geo_cancer_edgelist) <- c("GEO_Accession_Id","Cancer","Sample_Size")
geo_cancer_edgelist$GEO_Accession_Id <- as.character(as.vector(geo_cancer_edgelist$GEO_Accession_Id))
geo_cancer_edgelist$Cancer <- as.character(as.vector(geo_cancer_edgelist$Cancer))
geo_cancer_edgelist$Sample_Size <- as.numeric(as.vector(geo_cancer_edgelist$Sample_Size))

#Filter cancers with at least 100 samples 
geo_cancer_edgelist <- geo_cancer_edgelist[geo_cancer_edgelist$Sample_Size>100,]

final_expr_list <- list()
k <- 1
for (i in 1:nrow(geo_cancer_edgelist))
{
  geo_id <- geo_cancer_edgelist$GEO_Accession_Id[i]
  if (sum(names(expr_list1) %in% geo_id)>0)
  {
    id <- which(names(expr_list1)==geo_id)
    final_expr_list[[k]] <- expr_list1[[id]]
  }
  else if (sum(names(expr_list2) %in% geo_id)>0)
  {
    id <- which(names(expr_list2)==geo_id)
    final_expr_list[[k]] <- expr_list2[[id]]
  }
  k <- k + 1
}
names(final_expr_list) <- geo_cancer_edgelist$GEO_Accession_Id

load("../Data/Others/INV_genes.RData")

#Create the expression dataset and INV High and Low indices per cancer type
###############################################################################
for (i in 1:length(final_expr_list))
{
  #Get GEO Id and Cancer Id
  geo_id <- names(final_expr_list[i])
  cancer_type <- geo_cancer_edgelist[geo_cancer_edgelist$GEO_Accession_Id==geo_id,]$Cancer
  
  #Get sample ids of samples in dataset
  sample_ids <- as.character(rev_predictions_df[rev_predictions_df$Study==geo_id & rev_predictions_df$Cancer==cancer_type,]$ID)
  
  #Get the expression matrix from GEO accession
  if (i==3)
  {
    sample_names <- strsplit(colnames(final_expr_list[[i]]),"_")
    rev_sample_names <- NULL
    for (l in 1:length(sample_names))
    {
      rev_sample_names <- c(rev_sample_names,sample_names[[l]][1])
    }
    colnames(final_expr_list[[i]]) <- rev_sample_names
  }
  expr_matrix <- exprs(final_expr_list[[i]])[,sample_ids]
  
  #Map gene ids to gene names (HGNC symbols)
  gene_ids <- rownames(expr_matrix)
  
  if (i==6)
  {
    rev_gene_ids <- rep(0,length(gene_ids))
    for (l in 1:length(gene_ids))
    {
      gene_id <- unlist(strsplit(gene_ids[l]," "))[1]
      rev_gene_ids[l] <- as.numeric(gene_id)
    }
    query_out <- gconvert(query= rev_gene_ids, organism = "hsapiens", numeric_ns="ENTREZGENE_ACC",
                          target = "HGNC", mthreshold = Inf, filter_na = TRUE)
    gene_ids <- rev_gene_ids
  }else
  {
    query_out <- gconvert(query = gene_ids, organism = "hsapiens",
                          target="HGNC", mthreshold = Inf, filter_na = TRUE)
  }
  
  gene_names <- NULL
  for (j in 1:length(gene_ids))
  {
    if (gene_ids[j] %in% query_out$input)
    {
      gene_name <- query_out[query_out$input==gene_ids[j],]$name[1]
    }
    else{
      gene_name <- "NA"
    }
    gene_names <- c(gene_names,gene_name)
  }
  gene_indices <- which(gene_names!="NA" & gene_names!="nan")
  rev_expr_matrix <- expr_matrix[gene_indices,]
  rownames(rev_expr_matrix) <- gene_names[gene_indices]

  #Get the colMeans corresponding to INV genes
  inv_genes_in_matrix <- INV_genes[INV_genes %in% rownames(rev_expr_matrix)]
  inv_scores <- colMeans(rev_expr_matrix[inv_genes_in_matrix,])
  quantiles_inv <- quantile(inv_scores)
  quantile_low <- quantiles_inv[[2]]
  quantile_high <- quantiles_inv[[4]]
  inv_high_ids <- which(inv_scores>quantile_high)
  inv_low_ids <- which(inv_scores<quantile_low)

  #Get revised expression matrix with ordered ICR High and ICR Low
  D <- rev_expr_matrix[,c(inv_high_ids,inv_low_ids)]
  system(paste0("gunzip ../Results/",cancer_type,"/Adjacency_Matrix/",cancer_type,"_Full_Correlation_matrix.Rdata.gz"))
  load(paste0("../Results/",cancer_type,"/Adjacency_Matrix/",cancer_type,"_Full_Correlation_matrix.Rdata"))
  load(paste0("../Results/",cancer_type,"/Adjacency_Matrix/",cancer_type,"_Full_regulon_FGSEA.Rdata"))
  system(paste0("gzip ../Results/",cancer_type,"/Adjacency_Matrix/",cancer_type,"_Full_Correlation_matrix.Rdata"))
  
  #Revise the regulons with TFs present in expression matrix and targets present in expression matrix
  all_tfs <- names(regulon_sets)
  all_targets <- rownames(D)
  rev_tfs <- all_tfs[all_tfs %in% all_targets]
  rev_regulon_sets <- list()
  k <- 1
  for (tf in rev_tfs)
  {
    id <- which(all_tfs==tf)
    rev_targets <- regulon_sets[[id]][regulon_sets[[id]] %in% all_targets]
    rev_regulon_sets[[k]] <- rev_targets
    k <- k+1
  }
  names(rev_regulon_sets) <- rev_tfs
  rev_corr_matrix <- corr_matrix[,which(colnames(corr_matrix) %in% all_targets)]
  
  #Generate the activity matrix
  amat <- activity_mc(mexp = D, cormat = rev_corr_matrix, tflist = names(rev_regulon_sets), tau = 0.0)
  save(list=c('amat'),file=paste0("../Data/PRECOG/",cancer_type,"/",cancer_type,"_Full_Activity_matrix_FGSEA.Rdata"))
  inv_info <- c(rep("INV High",length(inv_high_ids)),rep("INV Low",length(inv_low_ids)))
  write.table(inv_info,paste0("../Data/PRECOG/",cancer_type,"/",cancer_type,"_INV_Info.csv"),row.names=T,col.names=F,quote=F,sep=",")
}

#Get the activity from multiple PRECOG datasets
########################################################################################################################
precog_cancers <- c("BLCA","BRCA","COAD","GBM","HNSC","LUAD","OV","SKCM")
inv_type <- c("INV Prognostic","INV Neutral","INV Prognostic","INV Neutral","INV Neutral","INV Prognostic","INV Prognostic","INV Neutral")
inputpath <- "../Data/PRECOG/"

cancer_info_vector <- NULL
inv_pheno_info_vector <- NULL
inv_type_info_vector <- NULL
sample_ids <- NULL

#Get total no of cancer samples
cancer_samples <- 0
for (i in 1:length(precog_cancers))
{
  cancer_type <- precog_cancers[i]
  table_indices <- read.table(paste0(inputpath,cancer_type,"/",cancer_type,"_INV_Info.csv"),header=FALSE,sep=",")
  cancer_samples <- cancer_samples+nrow(table_indices)
}

activity_matrix <- Matrix(0,nrow=length(all_mrs),ncol=cancer_samples)
rownames(activity_matrix) <- all_mrs
counter=0;

#Put values in the activity matrix
for (i in 1:length(precog_cancers))
{
  cancer_type <- precog_cancers[i]
  load(paste0(inputpath,cancer_type,"/",cancer_type,"_Full_Activity_matrix_FGSEA.Rdata"))
  table_indices <- read.table(paste0(inputpath,cancer_type,"/",cancer_type,"_INV_Info.csv"),header=FALSE,sep=",")
  table_indices$V2 <- as.character(as.vector(table_indices$V2))
  high_index_last <- sum(table_indices$V2=="INV High")
  low_index_last <- nrow(table_indices)
  
  amat_to_use <- amat  
  amat_to_use[amat_to_use>0] <- amat_to_use[amat_to_use>0]/max(amat_to_use)
  amat_to_use[amat_to_use<0] <- amat_to_use[amat_to_use<0]/abs(min(amat_to_use))
  samples <- dim(amat_to_use)[2]
  
  cancer_info_vector <- c(cancer_info_vector,rep(cancer_type,samples))
  inv_pheno_info_vector <- c(inv_pheno_info_vector,rep("INV High",sum(table_indices$V2=="INV High")),
                             rep("INV Low",sum(table_indices$V2=="INV Low")))
  inv_type_info_vector <- c(inv_type_info_vector,rep(inv_type[i],samples))
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
colnames(activity_matrix) <- sample_ids

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
activity_matrix <- as.matrix(activity_matrix)
#2nd column is INV High vs INV Low
hc.cols1 <- order(clab[,2],decreasing = T)
clab1 <- clab[hc.cols1,]
activity_matrix <- activity_matrix[,hc.cols1]
newhc.cols <- NULL
#3rd column is INV Enabled, INV Disabled and INV Neutral clusters
for (i in 1:length(unique(clab[,2])))
{
  color_id <- unique(clab[,2])[i]
  newhc.cols <- c(newhc.cols,length(newhc.cols)+order(clab1[clab1[,2]==color_id,3],decreasing = T))
}
activity_matrix <- activity_matrix[,newhc.cols]
clab2 <- clab1[newhc.cols,]
lower_val <- quantile(activity_matrix,probs=c(0.05))
upper_val <- quantile(activity_matrix,probs=c(0.95))
activity_matrix[activity_matrix<lower_val] <- lower_val
activity_matrix[activity_matrix>upper_val] <- upper_val

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
hc.cols.inv_high <- hclust(dist(t(activity_matrix[,all_inv_high_indices])),method='ward.D2')
hc.cols.inv_low <- hclust(dist(t(activity_matrix[,all_inv_low_indices])),method='ward.D2')
hc.cols <- c(hc.cols.inv_high$order,(length(all_inv_high_indices)+hc.cols.inv_low$order))
hc_col_list <- list(hc.cols.inv_high,hc.cols.inv_low)
hc_col <- .merge_hclust(hc_col_list)
#clab3 <- clab2[hc.cols,]
#Activity_Matrix <- Activity_Matrix[,hc.cols]

all_hot_mr_indices <- which(rlab[1,]==inv_colors[1])
all_cold_mr_indices <- which(rlab[1,]==inv_colors[2])
hc.rows.inv_high <- hclust(dist(activity_matrix[all_hot_mr_indices,]),method='ward.D2')
hc.rows.inv_low <- hclust(dist(activity_matrix[all_cold_mr_indices,]),method='ward.D2')
hc.rows <- c(hc.rows.inv_high$order,(length(all_hot_mr_indices)+hc.rows.inv_low$order))
hc_row_list <- list(hc.rows.inv_high,hc.rows.inv_low)
hc_row <- .merge_hclust(hc_row_list)

#Main Plotting Function for Figure 5B
pdf("../Results/Paper_Figures/Figure5/Heatmap_All_MRs_PRECOG_Figure_5B.pdf",height = 13, width=15, pointsize = 14)
par(bg="white")
par(fg="black",col.axis="black",col.main="black",col.lab="black", cex.main=2.0)
p1 <- heatmap.3(activity_matrix, scale="none", dendrogram="both", margins=c(4,20),
                Rowv=as.dendrogram(hc_row), Colv=as.dendrogram(hc_col), ColSideColors=clab2, labCol = FALSE,
                key=TRUE, keysize=1.5, labRow=all_mrs, cexRow=0.8, col=colorspace::diverge_hsv(100),RowSideColors=rlab,
                ColSideColorsSize=3, RowSideColorsSize=1, KeyValueName="Activity Value")
legend("topright",legend=c(unique(cancer_info_vector),"","",unique(inv_type_info_vector),"","",unique(inv_pheno_info_vector)),
       fill=c(unique(Cancers),"white","white",unique(INV_EorD),"white","white",unique(INV_Clusters)), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
dev.off()

ordered_mrs <- rev(rownames(activity_matrix)[p1$rowInd])
high_ids <- clab2[,2]=="red"
low_ids <- clab2[,2]=="blue"
activity_df <- perform_wilcox_test(activity_matrix[ordered_mrs,high_ids],activity_matrix[ordered_mrs,low_ids])
activity_df <- activity_df[,c(2:ncol(activity_df))]
activity_df$Mean1 <- round(activity_df$Mean1,3)
activity_df$Mean2 <- round(activity_df$Mean2,3)
activity_df$FC_Mean <- round(activity_df$FC_Mean,3)
activity_df$Pval <- signif(activity_df$Pval,3)
activity_df$Padjust <- signif(activity_df$Padjust,3)
activity_df <- activity_df[order(activity_df$Padjust,decreasing=T),]
write.table(activity_df,"../Results/Paper_Text/Diff_MRS_Activity_PRECOG_Table_S8.csv",row.names=T,col.names=T,quote=F,sep=",")

# #Make the heatmap for the set of deferentially active MRs identified for each cancer (INV-Low vs INV-High) on PRECOG datasets
# ###########################################################################################################################
# precog_cancers <- c("BLCA","BRCA","COAD","GBM","HNSC","LUAD","OV","SKCM")
# inputpath <- "../Data/PRECOG/"
# 
# for (i in 1:length(precog_cancers))
# {
#   cancer_type <- precog_cancers[i]
#   load(paste0(inputpath,cancer_type,"/",cancer_type,"_Full_Activity_matrix_FGSEA.Rdata"))
#   load(paste0(inputpath,cancer_type,"/",cancer_type,"_Full_TopMR_Info_FGSEA_BC_NES_1.Rdata"))
#   load(paste0(inputpath,cancer_type,"/",cancer_type,"_Full_Activity_matrix_FGSEA.Rdata"))
#   table_indices <- read.table(paste0(inputpath,cancer_type,"/",cancer_type,"_ICR_Info.csv"),header=FALSE,sep=",")
#   table_indices$V2 <- as.character(as.vector(table_indices$V2))
#   high_index_last <- sum(table_indices$V2=="ICR High")
#   low_index_last <- nrow(table_indices)
#   high_indices <- table_indices[table_indices$V2=="ICR High",]$V1
#   low_indices <- table_indices[table_indices$V2=="ICR Low",]$V1
#   
#   amat_to_use <- amat  
#   amat_to_use[amat_to_use>0] <- amat_to_use[amat_to_use>0]/max(amat_to_use)
#   amat_to_use[amat_to_use<0] <- amat_to_use[amat_to_use<0]/abs(min(amat_to_use))
#   samples <- dim(amat_to_use)[2]
#   
#   all_mrs_present <- topmr_info[topmr_info$pathway %in% rownames(amat_to_use),]$pathway
#   
#   #Get the activity matrix for these differential MRs per cancer with the phenotype information
#   new_mr_activity_matrix <- as.matrix(amat[all_mrs_present,c(high_indices,low_indices)])
#   rownames(new_mr_activity_matrix)
#   colnames(new_mr_activity_matrix) <- NULL
#   
#   colcol <- matrix(0,nrow=ncol(new_mr_activity_matrix),ncol=1)
#   high_ids <- c(1:length(high_indices));
#   low_ids <- c((length(high_indices)+1):length(colcol))
#   colcol[high_ids,1] <- "yellow"
#   colcol[low_ids,1] <- "green"
#   colnames(colcol) <- "ICR Phenotype"
#   hc_high.cols <- hclust(dist(t(new_mr_activity_matrix[,high_ids])),method='ward.D')
#   hc_low.cols <- hclust(dist(t(new_mr_activity_matrix[,low_ids])),method='ward.D')
#   hc.cols <- c(hc_high.cols$order,length(high_ids)+hc_low.cols$order);
#   hc.rows <- hclust(dist(new_mr_activity_matrix),method='ward.D')
#   
#   #Make the plot of activity matrix clustered by hot vs cold immune response
#   #=================================================================================================
#   pdf(paste0("../Results/",cancer_type,"/Images/",cancer_type,"_PRECOG_Activity_Matrix.pdf"),pointsize=9,height=9,width=15)
#   par(bg="white")
#   par(fg="black",col.axis="black",col.main="black",col.lab="black", cex.main=2.0, cex = 2.5)
#   new_mr_activity_matrix[new_mr_activity_matrix <= quantile(new_mr_activity_matrix, 0.05)] <- quantile(new_mr_activity_matrix, 0.05)
#   new_mr_activity_matrix[new_mr_activity_matrix >= quantile(new_mr_activity_matrix, 0.95)] <- quantile(new_mr_activity_matrix, 0.95)
#   heatmap.3(new_mr_activity_matrix[hc.rows$order,hc.cols], Rowv=TRUE, Colv=FALSE, col = bluered(100), scale="none", main= paste0("Master Regulator Activity for ",cancer_type," (PRECOG)"),
#             xlab = "TCGA Samples", ylab="Top MRs", dendrogram = "row", key = TRUE, density.info = "none", 
#             KeyValueName = "Activity Value", ColSideColors = colcol, ColSideColorsSize=2,
#             margins = c(6,6), useRaster = FALSE, cexRow = 0.6, cexCol = 2.0)
#   dev.off()
# }


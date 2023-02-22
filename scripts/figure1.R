####################################################################
###
### This Script clusters the RNASeq date from TCGA that have been
### normalized and filtered for specific tissue type
### and 1 sample per patient. This script includes a 
### log transformation of the data. Additionally, optimal kalinsky
### is calculated.
### 
### Input data:
### RNAseq data- 
### "Data/,"Cancer,"/",Cancer,"_gene_RNASeq_normalized_TP_filtered.Rdata"

### Output data are saved as Rdata file:
### "Data/",Cancer,"/",Cancer,"_INV_cluster_assignment_k2-6.Rdata"
### which includes: (1) table_cluster_assignment and (2) optimal.calinsky.
#####################################################################

rm(list=ls())

setwd("../")                                                                    # Setwd to location were output files have to be saved.

source("scripts/ipak.function.R")
source("scripts/get_functions.R")
source("scripts/High_Medium_Low_INV_classification.R")

required.packages = c("RCurl","httr", "rjson", "stringr", "HGNChelper","ggplot2")
required.bioconductor.packages = c("clue","ConsensusClusterPlus","survivalAnalysis","dplyr")
ipak(required.packages)
ibiopak(required.bioconductor.packages)
library(gplots)


#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

# Set Parameters
Cancers = c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA",
            "GBM","HNSC","KICH","KIRC","KIRP","LGG","LIHC","LUAD",
            "LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC",
            "SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

# Load INV genes data
load("Data/Others/INV_genes.RData") 
for (Cancer in Cancers)
{
  #Load RNA-Seq data
  out <- load(paste0("Data/",Cancer,"/",Cancer,"_gene_RNAseq_normalized_TP_filtered.Rdata"))
  D <- t(filtered.norm.RNAseqData)
  INV_subset_RNAseq = D[,colnames(D) %in% INV_genes]
  INV_subset_RNAseq_log2 = log(INV_subset_RNAseq +1, 2) 
  
  
  #Perform the consensus clustering
  ddist = dist(INV_subset_RNAseq_log2)
  
  ConsensusClusterObject <- ConsensusClusterPlus(ddist,                                                                
                                                 maxK = 6,                                                                              # set K
                                                 pItem = 0.8,
                                                 reps=5000,                                                                             # set repeats
                                                 title=paste0(Cancer,".INV.reps5000"),                    # Output filename
                                                 clusterAlg = "hc",                                                                     # clustering Algorithm : Hierarchiocal clustering
                                                 innerLinkage = "ward.D2",                                                              # for color coding the clusters use tmyPal = ...
                                                 finalLinkage = "complete",
                                                 plot = 'pdf',                                                                          # write resut to pdf (Alt:png)
                                                 writeTable = TRUE,
                                                 verbose = TRUE)
  outputfiles = list.files(paste0(Cancer,".INV.reps5000"), full.names = TRUE)
  class_files = outputfiles[grep("consensusClass", outputfiles)]
  
  N.files = length(class_files)
  table_cluster_assignment = data.frame(INVscore = rowMeans(INV_subset_RNAseq_log2))
  
  for (j in 1:N.files){
    file = paste0("./", class_files[j])
    consensus_class = read.csv(file = file,header=FALSE)
    group = paste0("Group_k",j+1)
    colnames(consensus_class) = c("PatientID", group)
    rownames(consensus_class) = consensus_class$PatientID
    consensus_class$PatientID = NULL
    table_cluster_assignment[,group] = consensus_class[,group][match(rownames(table_cluster_assignment), rownames(consensus_class))]
    
    transl_table_INV_cluster = aggregate(INVscore~get(group),data = table_cluster_assignment, FUN=mean)
    colnames(transl_table_INV_cluster) = c(group,"mean_INVscore")
    transl_table_INV_cluster = cbind(transl_table_INV_cluster[order(transl_table_INV_cluster$mean_INVscore),],INV_name=paste0("INV",c(1:(j+1))))
    
    INV_cluster = paste0("INV_cluster_k",j+1)
    table_cluster_assignment[,INV_cluster] = transl_table_INV_cluster$INV_name[match(table_cluster_assignment[,group],
                                                                                     transl_table_INV_cluster[,group])]
    
    
    
    
    #calinsky
    sHc <- hclust(ddist, method = "ward.D2")
    aCalinsky <- calinsky(sHc,gMax=10)
    pdf(file = paste0("Data/", Cancer, "/", Cancer, "_INV_cluster_assignment_k2-6.Calinsky.pdf"), width = 16, height = 6)
    plot(aCalinsky, type = "l", col = "grey", main = "Calinsky & Harabasz curve", xlab = "# of groups")
    text(1:length(aCalinsky), aCalinsky, paste(1:length(aCalinsky)))
    dev.off()
    optimal.calinsky = which(aCalinsky == max(aCalinsky[3:5]))
    
  }
  system(paste0("rm -rf ",Cancer,"*"))
  
  #save data
  outputname = paste0("Data/", Cancer, "/", Cancer, "_INV_cluster_assignment_k2-6_v2.Rdata")
  save(table_cluster_assignment,optimal.calinsky, file = outputname)
  table_cluster_assignment <- manual_annotation(outputname,Cancer)
  save(table_cluster_assignment,optimal.calinsky, file = outputname)
  
  #Make Heatmap highlighting the GSVA score for Invasiveness High vs Invasiveness Low vs Invasiveness Medium
  ##############################################################################
  temp_matrix <- scale(INV_subset_RNAseq_log2, center=T, scale=T)
  rev_INV_subset_RNAseq_log2 <- t(temp_matrix)
  inv_high_samples <- which(table_cluster_assignment$HML_cluster=="INV High")
  inv_low_samples <- which(table_cluster_assignment$HML_cluster=="INV Low")
  colcol <- matrix(0,nrow=ncol(rev_INV_subset_RNAseq_log2),ncol=1)
  colnames(colcol) <- "INV Pheno"
  
  #See if there are 2 or more clusters and prepare the heatmap accordingly
  ##############################################################################
  if (optimal.calinsky>2)
  {
    inv_medium_samples <- which(table_cluster_assignment$HML_cluster=="INV Medium")
    hc_high.cols <- hclust(dist(t(rev_INV_subset_RNAseq_log2[,inv_high_samples])),method='ward.D2')
    hc_low.cols <- hclust(dist(t(rev_INV_subset_RNAseq_log2[,inv_low_samples])),method='ward.D2')
    hc_medium.cols <- hclust(dist(t(rev_INV_subset_RNAseq_log2[,inv_medium_samples])),method='ward.D2')
    hc.cols <- c(hc_low.cols$order,length(inv_low_samples)+hc_medium.cols$order, length(inv_low_samples)+length(inv_medium_samples)+hc_high.cols$order);
    rev_INV_subset_RNAseq_log2 <- rev_INV_subset_RNAseq_log2[,c(inv_low_samples,inv_medium_samples, inv_high_samples)]
    colcol[hc_low.cols$order,1] <- "blue"
    colcol[length(hc_low.cols$order)+hc_medium.cols$order] <- "yellow"
    colcol[length(hc_low.cols$order)+length(hc_medium.cols$order)+hc_high.cols$order,1] <- "red"
  }
  hc.rows <- hclust(dist(rev_INV_subset_RNAseq_log2),method='ward.D2')
  
  #Save the Heatmap for each cancer 
  #############################################################################
  pdf(paste0("Results/",Cancer,"/",Cancer,"_Invasivenss_Expression.pdf"),height=4,width=4, pointsize=11, fonts = "sans")
  par(bg="white")
  par(fg="black",col.axis="black",col.main="black",col.lab="black", cex.main=1.0)
  heatmap.3(rev_INV_subset_RNAseq_log2[hc.rows$order,hc.cols], Rowv = NULL, Colv = NULL, col = bluered(110), scale="none", 
            main= paste0("Invasiveness Genes Expression for ",Cancer),
            xlab = "TCGA Samples", ylab = "", labRow = rownames(rev_INV_subset_RNAseq_log2)[hc.rows$order], labCol = FALSE, dendrogram = "none",
            key = TRUE, density.info = "none", KeyValueName = "Expression Value", ColSideColors = colcol, ColSideColorsSize=2, cexCol=2, 
            margins = c(5,5), useRaster = TRUE)
  dev.off()
  
}

#Perform the survival analysis w.r.t. PFI based on the invasiveness categories
###############################################################################
pb <- txtProgressBar(title = "Processing Cancers for Survival", min = 0, max = length(Cancers), width=300)
survival_list <- list()
for (i in 1:length(Cancers))
{
  Cancer <- Cancers[i] 
  
  #Get the curated survival information 
  clinical_df <- read.table(paste0("Data/survival_data/survival_",Cancer,"_survival.txt"),header=T,sep="\t")
  clinical_df <- clinical_df[!is.na(clinical_df$OS.time),]
  if (sum(is.na(clinical_df$Redaction)>0))
  {
    clinical_df <- clinical_df[is.na(clinical_df$Redaction),]  
  }else{
    clinical_df <- clinical_df[clinical_df$Redaction!="Redacted",]
  }
  
  #Load the table information with the clusters
  load(paste0("Data/",Cancer,"/",Cancer,"_INV_cluster_assignment_k2-6_v2.Rdata"))
  
  #Get the survival information w.r.t. the INV Score
  clinical_df$INVscore <- rep(NA, nrow(clinical_df))
  clinical_df$INV_Pheno <- rep(NA, nrow(clinical_df))
  clinical_df$TCGA_Match_Id <- rep(NA, nrow(clinical_df))
  rownames_table <- rownames(table_cluster_assignment)
  for (j in 1:nrow(clinical_df))
  {
    sample_id <- clinical_df[j,]$sample
    table_sample_id <- which(grepl(sample_id, rownames_table))
    if (length(table_sample_id)>0)
    {
      clinical_df[j,]$INVscore <- table_cluster_assignment[table_sample_id,]$INVscore
      clinical_df[j,]$INV_Pheno <- table_cluster_assignment[table_sample_id,]$HML_cluster
      clinical_df[j,]$TCGA_Match_Id <- rownames_table[table_sample_id]
    }
  }
  clinical_df <- clinical_df[!is.na(clinical_df$INVscore),]
  clinical_df  <- clinical_df[!is.na(clinical_df$TCGA_Match_Id),]
  clinical_df <- clinical_df[clinical_df$OS.time<=3650,]
  clinical_df <- clinical_df[clinical_df$INV_Pheno!="INV Medium",]
  
  ##Find the optimal cutpoint based on the INV Score
  #res.cut <- surv_cutpoint(data=clinical_df, time="PFI.time", event="PFI", variables = "INVscore", minprop=0.1, progressbar = T)
  #res.cat <- surv_categorize(res.cut, variables = "INVscore", labels=c("INV Low","INV High"))
  #clinical_df$INV_Pheno <- res.cat$INVscore
  
  #Make the survival plot based on High vs Low categories
  print(Cancer)
  fit <- analyze_survival(data=clinical_df,time_status = vars(OS.time, OS), by = INV_Pheno)
  print(fit)
  survival_list[[i]] <- fit
  g <- kaplan_meier_plot(fit, 
                         break.time.by="breakByYear",
                         xscale = 365.25,
                         risk.table = F,
                         pval = F,
                         hazard.ratio=F,
                         xlab = "Time (in Years)",
                         palette=c("red","blue"),
                         ggtheme=ggplot2::theme_bw(base_size=10))
  g$plot <- g$plot+theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         text = element_text(size=10, colour="black" ),
                         axis.text.x = element_text(size=10,  colour = "black"),
                         axis.text.y = element_text(size=10,  colour = "black" ),
                         axis.title.x = element_text(size=10,  colour = "black"),
                         axis.title.y = element_text(size=10,  colour = "black"),
                         plot.title = element_text(size = 10,  colour="black"),
                         axis.text = element_text(size=10,  colour="black"))
  ggsave(filename=paste0("Results/",Cancer,"/",Cancer,"_Survival_Categorization.pdf"), plot = print(g$plot, newpage = FALSE), device = cairo_pdf(pointsize = 10), height=2, width=3.5, units="in", dpi=300)
  dev.off()
  write.table(clinical_df, file=paste0("Results/",Cancer,"/",Cancer,"_survival_info.csv"),row.names=F, col.names=T, sep="\t", quote=F)
  
  setTxtProgressBar(pb, i, label=paste0("Done for Cancer: ",Cancer))
  
}
names(survival_list) <- Cancers

#Make the forest plot data frame
########################################################################
library(dplyr)
library(tcltk)
library(forestplot)

mean <- NULL
UI <- NULL
LI <- NULL
cancer_info <- NULL
pval_info <- NULL
hazards_info <- NULL
invhigh_info <- NULL
invlow_info <- NULL
for (i in 1:length(survival_list)){
  cancer <- names(survival_list)[i]
  fit_info <- survival_data_frames(survival_list[[i]])
  temp_mean <- as.numeric(fit_info$hazardRatios[2,"HR"])
  temp_ui <- as.numeric(fit_info$hazardRatios[2,"Upper.CI"])
  temp_li <- as.numeric(fit_info$hazardRatios[2,"Lower.CI"])
  if (!is.infinite(temp_ui) & !is.infinite(temp_li))
  {
    temp_pval <- as.character(fit_info$hazardRatios[2,"p"])
    temp_invhigh <- as.character(fit_info$strataMetadata$records[1])
    temp_invlow <- as.character(fit_info$strataMetadata$records[2])
    if (as.numeric(temp_invhigh)+as.numeric(temp_invlow)>=0)
    {
      if (temp_pval=="<0.001")
      {
        temp_pval <- paste0(temp_pval," ***")
      }
      else{
        temp_pval <- as.numeric(temp_pval)
        if (temp_pval<=0.1 & temp_pval>0.05)
        {
          temp_pval <- paste0(round(temp_pval,2)," *")
        }
        else if (temp_pval<=0.05 & temp_pval>=0.001)
        {
          temp_pval <- paste0(temp_pval," **")
        }
        else if (temp_pval>0.1)
        {
          temp_pval <- paste0(temp_pval)
        }
      }
      mean <- c(mean, temp_mean)
      UI <- c(UI, temp_ui)
      LI <- c(LI, temp_li)
      pval_info <- c(pval_info, temp_pval)
      cancer_info <- c(cancer_info, cancer)
      invhigh_info <- c(invhigh_info, temp_invhigh)
      invlow_info <- c(invlow_info, temp_invlow)
      hazards_info <- c(hazards_info, paste0(temp_mean,"(",temp_li,"-",temp_ui,")"))
    }
  }
}
base_data <- tibble(Cancer = cancer_info,
                    N.INV.High = invhigh_info,
                    N.INV.Low = invlow_info,
                    P.value = pval_info,
                    HR = hazards_info,
                    mean = mean,
                    upper = UI,
                    lower = LI,
)
base_data <- base_data[order(mean, decreasing=F),]
header <- tibble(Cancer = "Cancer",
                 N.INV.High = "N1",
                 N.INV.Low = "N2",
                 P.value = "P-value",
                 HR = "HR(95%CI)"
)
output_df <- bind_rows(header,base_data)

#Make the plot for all cancers
fn_uni_full <- local({
  i = 0
  b_clrs = c(rep("blue",3),rep("red",29))
  
  function(..., clr.line, clr.marker){
    i <<- i + 1
    fpDrawNormalCI(..., clr.line = b_clrs[i], clr.marker = b_clrs[i])
  }
})

output_df %>% 
  forestplot(labeltext = c(Cancer, N.INV.High, N.INV.Low, P.value, HR), 
             fn.ci_norm = fn_uni_full,
             graph.pos = 6,
             boxsize = 0.25,
             xlog = TRUE,
             graphwidth = unit(4,"cm"),
             colgap = unit(10,"mm"),
             #col = fpColors(box = "blue",line = "blue"),
             lwd.xaxis=1,
             hrzl_lines = list("2" = gpar(lty = 2), 
                               "5" = gpar(lwd = 1, columns = 1:6, col = "#000044")),
             lineheight = unit(5,"mm"),
             txt_gp = fpTxtGp(ticks=gpar(fontfamily="sans",fontsize=20),
                              label=gpar(fontfamily="sans", fontsize=10),
                              title=gpar(fontfamily="sans", fontsize=10, fontface="bold"))) -> p_full
pdf("Results/PanCancer/All_Cancers_of_Interest_Invasivenss_Forest_Plot_Supp.pdf", width = 8, height=8, pointsize = 11)
p_full
dev.off()

#Make the Boxplot highlight the PANoptosis score for all cancers
################################################################################
boxplot_df <- NULL
diff_mean_low_high <- NULL
for (Cancer in Cancers)
{
  load(paste0("Data/",Cancer,"/",Cancer,"_INV_cluster_assignment_k2-6_v2.Rdata"))
  mean_low <- mean(table_cluster_assignment[table_cluster_assignment$HML_cluster=="INV Low",]$INVscore)
  mean_high <- mean(table_cluster_assignment[table_cluster_assignment$HML_cluster=="INV High",]$INVscore)
  diff_mean_low_high <- c(diff_mean_low_high,mean_high)
  temp <- cbind(rep(Cancer,nrow(table_cluster_assignment)),table_cluster_assignment$INVscore,table_cluster_assignment$HML_cluster)
  boxplot_df <- rbind(boxplot_df, temp)
}
boxplot_df <- as.data.frame(boxplot_df)
colnames(boxplot_df) <- c("Cancer","INVscore","INV_Category")
boxplot_df$Cancer <- as.character(as.vector(boxplot_df$Cancer))
boxplot_df$INVscore <- as.numeric(as.vector(boxplot_df$INVscore))
boxplot_df$INV_Category <- as.character(as.vector(boxplot_df$INV_Category))
boxplot_df$INV_Category <- factor(boxplot_df$INV_Category, levels=c("INV Low","INV Medium","INV High"))
names(diff_mean_low_high) <- Cancers
diff_mean_low_high <- sort(diff_mean_low_high, decreasing = F)
boxplot_df$Cancer <- factor(boxplot_df$Cancer, levels = names(diff_mean_low_high))


g_boxplot <- ggplot(data=boxplot_df, aes(x=Cancer, y=INVscore, fill=INV_Category))+
  geom_boxplot(position=position_dodge(1)) + scale_fill_manual(name="", values=c("blue","yellow","red"))+
  theme_bw() + xlab("Cancers") + ylab("Invasiveness Score")+ylim(c(4,16))+
  theme(plot.background=element_rect(fill = "white"),
        panel.background = element_rect(fill = 'white'),legend.background = element_rect(fill = "white", color = NA),
        legend.key = element_rect(color = "gray", fill = "white"), legend.title = element_text(color = "black"), legend.text = element_text(size=10, family="sans", color = "black"),
        axis.text.x = element_text(size=10, angle=90, hjust = 0.5, vjust=0.5, family="sans", color="black"),
        axis.text.y = element_text(size=10, family="sans",color="black"),
        axis.title.x = element_text(size=10, family="sans", color="black"),
        axis.title.y = element_text(size=10, family="sans", color="black"))
ggsave(filename = "Results/PanCancer/All_Cancers_of_Interest_Boxplot_Invasiveness_Score.pdf", plot = g_boxplot, device = pdf(fonts="sans"), width = 9, height=4, units="in", dpi=300)
dev.off()



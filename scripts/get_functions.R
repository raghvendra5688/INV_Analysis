library(devtools)
#install_github("miccec/yaGST")
require("yaGST")
require("parallel")

#Functions

#============================================================================================
loading_data <- function(filename,M){
  z <- readRDS(file=paste0("../Data/",filename,"/TCGA-",filename,"_normcounts.rda"))
  
  #Get the gene names and associate it with appropriate variables
  #====================================================================================
  load("../Data/Others/TF_Target_Info.Rdata")
  tfs <- tf_gene_names;
  targets <- target_gene_names;
  rownames(z) <- target_gene_names;

  #Get TFs and targets
  all_genes <- rownames(z);
  m_targets <- colnames(M);
  m_tfs <- rownames(M);
  temp_target_genes <- intersect(all_genes,m_targets);
  tf_genes <- intersect(all_genes,m_tfs);
  rem_target_genes <- setdiff(temp_target_genes,tf_genes);
  target_genes <- c(tf_genes,rem_target_genes);
  
  #Get the modified gene expression matrix
  genes_from_expression_matrix <- union(target_genes,tf_genes);
  modified_exp_matrix <- t(z[genes_from_expression_matrix,]);
  rm(z);
  gc();
  N <- nrow(modified_exp_matrix);
  d <- ncol(modified_exp_matrix);
  
  #Get final mechanistic network
  g_M <- M[tf_genes,target_genes];
  
  #Make the perturbation matrix
  K <- matrix(0,nrow=N,ncol=d);
  colnames(K) <- colnames(modified_exp_matrix);
  
  return(list(modified_exp_matrix,K,g_M,tf_genes,target_genes));
}

#===========================================================================================
get_high_low_indices <- function(table_cluster_assignment, D){
  
  #Function to get all the ICR High and ICR low samples (i.e. phenotype information from TCGA)
  tcga_sample_ids <- rownames(table_cluster_assignment)
  unmatched_tcga_sample_ids <- colnames(D);
  filter_samples <- which(unmatched_tcga_sample_ids %in% tcga_sample_ids);
  D <- D[,filter_samples];
  order_indices <- NULL
  for (i in 1:length(unmatched_tcga_sample_ids))
  {
    order_indices <- c(order_indices,which(tcga_sample_ids==unmatched_tcga_sample_ids[i]));
  }
  table_cluster_assignment <- table_cluster_assignment[order_indices,];
  high_indices <- which(table_cluster_assignment$HML_cluster=="INV High")
  low_indices <- which(table_cluster_assignment$HML_cluster=="INV Low")
  return(list(table_cluster_assignment,high_indices,low_indices,D))
}

#===========================================================================================
parcor <- function(Mat) { 
  
  #Calculate correlation matrix in parallel
  require(parallel) 
  nc <- detectCores() 
  m <- dim(Mat)[1] 
  n <- dim(Mat)[2] 
  rL <- split(t(Mat),1:n) 
  res <- mclapply(rL, function(x){ mm <- matrix(unlist(x), ncol = m, byrow = TRUE) 
                  cor(t(mm), Mat)
                  },
                  mc.cores = nc ) 
  out <- matrix(unlist(res), nrow=n, byrow = T) 
  return(out)
}


#==========================================================================================
get_regulons <- function(net,corr_matrix,minsize=20){
  regulons=list()
  tnames=c()
  
  for(tfi in rownames(net)){
    tgs <- which(net[tfi,]!=0)
    pos<- which(corr_matrix[tfi,tgs]>0)
    neg<- which(corr_matrix[tfi,tgs]<0)
    r=colnames(net[tfi,net[tfi,]!=0,drop=F])
    if(length(r)>minsize & length(pos)>1 & length(neg)>1) {
      tnames=c(tnames,tfi)
      regulons[[length(regulons)+1]] = r
    }
  }
  names(regulons)=tnames
  return(regulons)
}

#====================================================================================================
get_viper_network <- function(net, minsize)
{
  #Prepare the 3 column sparse gene regulatory network for Viper
  net_sparse <- as(as.matrix(net), "dgCMatrix")
  net_df <- as.data.frame(summary(net_sparse))
  net_df$TF <- rownames(net_sparse)[net_df$i]
  net_df$Target <- colnames(net_sparse)[net_df$j]
  net_df <- net_df[,-c(1,2)]
  net_df <- net_df[,c(2,3,1)]
  colnames(net_df) <- c("TF","Target","MI")
  revised_net_df <- NULL
  
  #Filtering based on Minimum no of targets regulated by a TF
  for (tf in tfs)
  {
    temp_edges <- net_df[net_df$TF==tf,]
    if (nrow(temp_edges)>minsize)
    {
      revised_net_df <- rbind(revised_net_df,temp_edges)
    }
  }
  revised_net_df <- as.data.frame(revised_net_df)
  colnames(revised_net_df) <- c("TFs","Targets","Weight")
  revised_net_df$TFs <- as.character(as.vector(revised_net_df$TFs))
  revised_net_df$Targets <- as.character(as.vector(revised_net_df$Targets))
  revised_net_df$Weight <- as.numeric(as.vector(revised_net_df$Weight))
  return(revised_net_df)
}


#==========================================================================================
#Perform wilcox test on two expression matrices for two cases using identified tfs
perform_wilcox_test <- function(A,B,exact=FALSE)
{
  genes <- rownames(A);
  wilcox_test_info <- NULL;
  for (i in 1:length(genes))
  {
    temp <- wilcox.test(A[i,],B[i,],exact=FALSE)
    statistic <- temp$statistic
    p_value <- temp$p.value
    if (is.nan(p_value)) { p_value = 1};
    m1 <- mean(A[i,]);
    m2 <- mean(B[i,]);
    FC_mean <- m1-m2;
    wilcox_test_vector <- c(statistic,p_value,m1,m2,FC_mean);
    wilcox_test_info <- rbind(wilcox_test_info,wilcox_test_vector);
  }
  wilcox_test_info <- as.data.frame(wilcox_test_info);
  colnames(wilcox_test_info) <- c("Statistic","Pval","Mean1","Mean2","FC_Mean")
  wilcox_test_info$Statistic <- as.numeric(as.vector(wilcox_test_info$Statistic))
  wilcox_test_info$Pval <- as.numeric(as.vector(wilcox_test_info$Pval))
  wilcox_test_info$Mean1 <- as.numeric(as.vector(wilcox_test_info$Mean1))
  wilcox_test_info$Mean2 <- as.numeric(as.vector(wilcox_test_info$Mean2))
  wilcox_test_info$FC_Mean <- as.numeric(as.vector(wilcox_test_info$FC_Mean))
  wilcox_test_info$Padjust <- p.adjust(wilcox_test_info$Pval,method="fdr")
  rownames(wilcox_test_info) <- genes;
  return(wilcox_test_info)
}

#==========================================================================================
#Get the activity matrix
activity_mc <- function(mexp,cormat,tflist=NULL,tau=0.6) {

  mexp.s = mexp
  for(i in 1:nrow(mexp.s)){
    #mexp.s[i,] = (mexp.s[i,] - mean(mexp[i,]))/sd(mexp.s[i,])
    mexp.s[i,] = (mexp.s[i,]-mean(mexp[i,]));
  }
  if (is.null(tflist)) {
    tflist=rownames(cormat)
  }
  
  actmat = mexp.s[tflist,]
  actmat[1:length(actmat)]=0
  pb = txtProgressBar(min=1,max=length(tflist),style=3)
  i=1
  for(tfi in tflist) {
    postrg = names(cormat[tfi,cormat[tfi,] > tau])
    negtrg= names(cormat[tfi,cormat[tfi,] < -tau])
    if (length(postrg)>1 & length(negtrg)>1) {
      apos = apply(mexp.s[postrg,,drop=F],2,sum)/length(postrg)
      aneg = apply(mexp.s[negtrg,,drop=F],2,sum)/length(negtrg)
      actmat[tfi,] =  apos - aneg
    } 
    else if (length(postrg)>1 & length(negtrg)<=1)
    {
      actmat[tfi,] = 1;
    }
    else if (length(postrg)<=1 & length(negtrg)>1)
    {
      actmat[tfi,] = -1;
    }
    else
    {
      actmat[tfi,] = 0;
    }
    
    setTxtProgressBar(pb, i)
    i=i+1

  }
  return(actmat)
}

#==========================================================================================
#Get the activity matrix considering tumor purity
activity_mc_with_purity <- function(mexp,cormat,tflist=NULL,tau=0.6,purity) {
  
  mexp.s = mexp * purity
  for(i in 1:nrow(mexp.s)){
    #mexp.s[i,] = (mexp.s[i,] - mean(mexp[i,]))/sd(mexp.s[i,])
    mexp.s[i,] = (mexp.s[i,]-mean(mexp[i,]));
  }
  if (is.null(tflist)) {
    tflist=rownames(cormat)
  }
  actmat = mexp.s[tflist,]
  actmat[1:length(actmat)]=0
  pb = txtProgressBar(min=1,max=length(tflist),style=3)
  i=1
  for(tfi in tflist) {
    postrg = names(cormat[tfi,cormat[tfi,] > tau])
    negtrg= names(cormat[tfi,cormat[tfi,] < -tau])
    if (length(postrg)>1 & length(negtrg)>1) {
      apos = apply(mexp.s[postrg,,drop=F],2,sum)/length(postrg)
      aneg = apply(mexp.s[negtrg,,drop=F],2,sum)/length(negtrg)
      actmat[tfi,] =  apos - aneg
    } 
    else if (length(postrg)>1 & length(negtrg)<=1)
    {
      actmat[tfi,] = 1;
    }
    else if (length(postrg)<=1 & length(negtrg)>1)
    {
      actmat[tfi,] = -1;
    }
    else
    {
      actmat[tfi,] = 0;
    }
    
    setTxtProgressBar(pb, i)
    i=i+1
    
  }
  return(actmat)
}


#===========================================================================================
#Another function to create regulons given adjacency matrix of regulatory network, correlation matrix between TF-target and minsize
createRegulons <- function(net,ccor,minsize=20)
{
  regulon <- vector("list",nrow(net))
  for(i in 1:nrow(net)){
    tgs <- which(net[i,]!=0)
    pos <- which(ccor[i,tgs]>0)
    neg <- which(ccor[i,tgs]<0)
    if((length(pos)>minsize | length(neg)> minsize) & (length(pos)>1 & length(neg)>1)){
      regulon[[i]] = list(pos=names(pos),neg=names(neg))
    }
    else{
      regulon[[i]] = NA
    }
  }
  names(regulon)<- rownames(net)
  regulon <- regulon[!is.na(regulon)]
  return(regulon)
}

#============================================================================================
#Infer gene regulatory network using ARACNE type approach
massiveGST <- function(rrnk, GGO, alternative = "greater", keepDetails = FALSE, 
                       writeXLS = FALSE, fName = "massiveGST.xls") {
  size <- length(GGO)
  GO <-  intersect(GGO, names(rrnk))
  actualSize <- length(GO)
  result <- data.frame(collection=c("Set"), size, actualSize)
  rnk <- rank(rrnk)
  sumOfRanks <-  sum(rnk[GO])
  n <- length(rnk)
  nx <- actualSize #ny
  ny <- n - nx #nx
  U_stat <- nx * ny + ny * (ny + 1)/2 + sumOfRanks - n * (n + 1)/2
  pod <- U_stat/nx/ny 
  odd <- pod/(1-pod)
  log2_odd <- log2(odd)
  zValue <- U_stat - nx * ny/2
  sigma <- sqrt(nx * ny * (nx + ny + 1)/12)
  correction <- switch(alternative, two.sided = sign(zValue) * 0.5, greater = 0.5, less = -0.5)
  zValue <- (zValue - correction)/sigma
  pValue <- switch(alternative, less = 1 - pnorm(zValue), greater = pnorm(-zValue), two.sided = 2 * pnorm(-abs(zValue)))
  qValue <- p.adjust(pValue, method = "BH") # fdr
  result <- data.frame(result, sumOfRanks, U_stat, pod, odd, log2_odd, zValue, pValue, qValue)
  rowsToRemove <- which(result[, "actualSize"] == 0)
  if(length(rowsToRemove) > 0) result <- result[-rowsToRemove,]
  if(alternative == "two.sided") {
    abs_log2_odd <- abs(result[, "log2_odd"])
    result <- cbind(result, abs_log2_odd)
    positive <- which(result[, "log2_odd"] > 0)
    order_p <- (rank(result[positive, "log2_odd"]) + rank(result[positive, "actualSize"]) + rank(-log10(result[positive, "pValue"])))/3
    names(order_p) <- rownames(result)[positive]
    negative <- which(result[, "log2_odd"] <= 0)
    order_n <- (rank(-result[negative, "log2_odd"]) + rank(result[negative, "actualSize"]) + rank(- log10(result[negative, "pValue"])))/3
    names(order_n) <- rownames(result)[negative]
    ordering <- c(order_p, -order_n)[rownames(result)]
    result <- cbind(result[, 1:12], ordering)
  } else {
    if(alternative == "greater") {
      ordering <- (rank(result[, "log2_odd"]) + rank(result[, "actualSize"]) + rank(-log10(result[, "pValue"])))/3
      names(ordering) <- rownames(result)[rownames(result)]
      result <- cbind(result, ordering)
    } else {
      ordering <- (rank(-result[, "log2_odd"]) + rank(result[, "actualSize"]) + rank(-log10(result[, "pValue"])))/3
      names(ordering) <- rownames(result)[rownames(result)]
      result <- cbind(result, ordering)
    }
  }
  
  ordering <- order(result[, "ordering"], decreasing = TRUE)
  result <- result[ordering,]
  ##########
  if(!keepDetails) {
    colsToRemove <- c("sumOfRanks", "U_stat", "zValue")
    colsToRemove <- which(colnames(result) %in% colsToRemove)
    result <- result[, -colsToRemove]
  }
  
  if(writeXLS) {
    gstTable <- as.data.frame(result)
    require(WriteXLS)
    WriteXLS("gstTable", ExcelFileName = fName, row.names = TRUE)
  }
  invisible(result)
}

massiveExtGst <- function (rankedList, geneSetUp, geneSetDown, minLenGeneSet = 15) 
{
  
  doubleRankedList <- c(rankedList, -rankedList)
  flag <- c(rep(TRUE, length(rankedList)), rep(FALSE, length(rankedList)))
  oorder <- order(doubleRankedList, decreasing = TRUE)
  doubleRankedList <- doubleRankedList[oorder]
  flag <- flag[oorder]
  hits <- (names(doubleRankedList) %in% geneSetUp) & flag
  hits <- hits | ((names(doubleRankedList) %in% geneSetDown) & 
                    !flag)
  names(doubleRankedList)[which(!hits)] <- paste0(names(doubleRankedList)[which(!hits)], 
                                                  1:sum(!hits))
  geneSet <- c(geneSetUp, geneSetDown)
  ans <- massiveGST(doubleRankedList, geneSet,alternative = "two.sided")
  invisible(ans)
}


#=============================================================================================
#Compute activity based on network inferred using the massiveGST
computeActivity <- function(net,E,regulons,ncores=16){
  M<- apply(E,1,mean)
  E <- E-M
  tf<- intersect(rownames(net),rownames(E))
  targets  <- intersect(colnames(net),rownames(E))
  net <- net[tf,targets]
  E <- E[targets,]
  
  aMat <- matrix(0,nrow=length(regulons),ncol=ncol(E))
  
  for(i in 1:ncol(E)){
    rankedList <- sort(E[,i], decreasing=T)
    #lans <- mclapply(regulon,function(x) mwwExtGST(rankedList = rankedList,geneSetUp = x$pos,
    #                                             geneSetDown = x$neg),mc.cores = ncores)
    lans <- mclapply(regulons,function(x) massiveExtGst(rankedList = rankedList,geneSetUp = x$pos,
                                                       geneSetDown = x$neg),mc.cores = ncores)
    
    nes <- unlist(lapply(lans,function(x) x$pod ))
    activity <- log2(nes/(1-nes))
    aMat[,i] <- activity
  }
  return(aMat)
}

######################
## Stefano's functions
##
####################

required.packages <- c("clue")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]

if(length(missing.packages)) install.packages(missing.packages)


## Quantile Normalization
quantileNormalization <- function (wd, distribution) 
{
  n <- nrow(wd)
  m <- ncol(wd)
  if (!missing(distribution)) 
    if (m != length(distribution)) 
      stop("The reference distribution has length different from the col dimension of the data matrix.")
  else distribution <- sort(distribution)
  o <- matrix(0, n, m)
  for (i in 1:n) o[i, ] <- order(wd[i, ])
  j <- 1
  tmp <- rep(0, n)
  while (j <= m) {
    if (missing(distribution)) {
      for (i in 1:n) tmp[i] <- wd[i, o[i, j]]
      value <- mean(tmp)
    }
    else value <- distribution[j]
    for (i in 1:n) wd[i, o[i, j]] <- value
    j <- j + 1
  }
  return(wd)
}


## Calinsky's function to determine optimal k value
calinsky <- function (hhc, ddist = NULL, gMax = round(1 + 3.3 * log(length(hhc$order), 
                                                                    10))) 
{
  ans <- rep(0, gMax)
  uDist <- as.cl_ultrametric(hhc)
  uDist <- as.matrix(uDist)
  totalMedoid <- which.min(apply(uDist, 1, sum))
  totalSum <- sum(uDist[totalMedoid, ])
  n <- length(hhc$order)
  for (g in 2:gMax) {
    cclust <- cutree(hhc, k = g)
    groupMedoids <- rep(0, g)
    withinSum <- 0
    for (k in 1:g) {
      if (sum(cclust == k) == 1) 
        next
      tmp <- which.min(apply(uDist[cclust == k, cclust == 
                                     k], 1, sum))
      groupMedoids[k] <- which(names(cclust) == names(tmp))
      withinSum <- withinSum + sum(uDist[groupMedoids[k], 
                                         cclust == k])
    }
    betweenSum <- totalSum - withinSum
    betweenSum <- betweenSum/(g - 1)
    withinSum <- withinSum/(n - g)
    ans[g] <- betweenSum/withinSum
  }
  class(ans) <- "calinsky"
  return(ans)
}

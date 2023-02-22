library(corto)
library(doMC)
library(foreach)
library(igraph)
registerDoMC(20)

setwd('../scripts/')
source('get_functions.R')
source('gene-reverse-network.R')

#Load mechanistic network
load('../Data/Others/me_net_full.Rdata')

#Do the GRN construction and saving
#========================================================================================
filename <- "MESO"
outputpath <- paste0("../Results/ARACNE/",filename)
sample_name <- paste0(filename,"_Full_");
out <- loading_data(filename,M)
D <- t(log2(out[[1]]+1))
tfs <- rownames(M)

#Get network using ARACNE-AP algorithm with number of bootstraps to 1000 and perform DPI too
net <- corto(D, tfs, nbootstraps = 100, nthreads = 20, p=1e-3, verbose=TRUE)

V_final <- matrix(0,nrow=length(net),ncol=nrow(D))
rownames(V_final) <- names(net)
colnames(V_final) <- rownames(D)
for (i in 1:length(net))
{
  tf <- names(net)[i]
  targets <- names(net[[i]]$tfmode)
  V_final[tf,targets] <- net[[i]]$likelihood
}
print(sum(V_final>0))
write.table(V_final,file=paste0(outputpath,"/Adjacency_Matrix/Final_Adjacency_Matrix_v2.csv"),row.names = T, col.names = T)


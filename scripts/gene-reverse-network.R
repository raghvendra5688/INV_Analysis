aracne2 = function(mexp,from,to,nperm=1000) {
  require(parmigene)
  pb = txtProgressBar(min=1,max=nperm,style=3)
  mi = knnmi.cross(mexp[from,],mexp[to,])
  mi.perm = array(0,dim=c(nrow(mi),ncol(mi)))
  rates = array(0,dim=c(nrow(mi),ncol(mi)))
  for (i in 1:nperm) {
    mi.tmp = knnmi.cross(mexp[from,],mexp[to,sample(1:ncol(mexp))])
    mi.perm = mi.perm + (mi < mi.tmp)
    rates = rates + abs(mi.tmp)
    setTxtProgressBar(pb, i)
  }
  
  pval = mi.perm/nperm # higher tail
  rates = nperm/rates
  
  mi0 = mi[pval==0]
  ra0 = rates[pval==0]
  pval[pval==0] = sapply(1:length(mi0),function(i) pexp(mi0[i],rate=ra0[i],lower.tail = F))
  
  pval.fdr = p.adjust(pval,method = "fdr")
  mi[pval.fdr<=0.05] = 0
 
  return(list(MI=mi,PVAL=pval,NPERM=nperm))
}

dpi2 = function(mi,from,to) {
  require(parmigene)
  commong = intersect(from,to)
  onlyfrom = setdiff(from,to)
  onlyto = setdiff(to,from)
  extn = length(commong)+length(onlyfrom)+length(onlyto)
  extmi = matrix(0,nrow=extn,ncol=extn)
  rownames(extmi) = c(commong,onlyto,onlyfrom)
  colnames(extmi) = c(commong,onlyto,onlyfrom)
  extmi[commong,commong] = mi[commong,commong,drop=FALSE]
  extmi[onlyto,commong] = t(mi[commong,onlyto,drop=FALSE])
  extmi[onlyfrom,commong] = mi[onlyfrom,commong,drop=FALSE]
  extmi[commong,onlyfrom] = t(mi[onlyfrom,commong,drop=FALSE])
  extmi[onlyto,onlyfrom] = t(mi[onlyfrom,onlyto,drop=FALSE])
  extmi[commong,onlyto] = mi[commong,onlyto,drop=FALSE]
  extmi[onlyfrom,onlyto] = mi[onlyfrom,onlyto,drop=FALSE]
  extmi = aracne.a(extmi)
  return(extmi[from,to])
}


cross.cor <- function(x, y, verbose = TRUE, ncore="all", met="pearson", ...){
  require(doParallel)
  if(ncore=="all"){
    ncore = parallel:::detectCores()
    doParallel:::registerDoParallel(cores=ncore)
  } else{
    doParallel:::registerDoParallel(cores=ncore)
  }
  N <- nrow(x)
  M <- nrow(y)
  
  ntasks <- ncore
  TM<-round(sqrt(ntasks*M/N))
  TN<-ntasks %/% TM
  if (TM==0) TM=1
  if (TN==0) TN=1
  
  Nsize <- N %/% TN
  Msize <- M %/% TM
  corMAT<-foreach(i = 1:TN, .combine='rbind') %:% 
    foreach(j = 1:TM, .combine='cbind') %dopar% {
      s1<-(i-1)*Nsize+1
      e1<-s1+Nsize-1
      s2<-(j-1)*Msize+1
      e2<-s2+Msize-1     
      if(i==TN) {
        e1<-N
      }
      if(j==TM) {
        e2<-M
      }
      #cat(s1,e1,s2,e2,"\n")
      cor(t(x[s1:e1,]), t(y[s2:e2,]), method = met, use="pairwise.complete.obs", ...)
    }
  gc()
  corMAT[is.na(corMAT)]<- 0;
  return(corMAT)
}


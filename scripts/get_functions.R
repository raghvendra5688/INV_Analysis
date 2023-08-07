library(devtools)
#install_github("miccec/yaGST")
require("yaGST")
require("parallel")

#Functions
heatmap.3 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.3,
                      cexRow = 0.2 + 1/log10(nr),
                      cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      ColSideColorsSize = 1,
                      RowSideColorsSize = 1,
                      KeyValueName="Value",...){
  
  invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid)))
    else if (is.vector(x))
      return(all(is.na(x)))
    else return(FALSE)
  }
  
  x <- as.matrix(x)
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(dim(Rowv)) || is.na(dim(Rowv)))
    Rowv <- FALSE
  if (is.null(dim(Colv)) || is.na(dim(Colv)))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    
    if (!missing(ColSideColors)) {
      #if (!is.matrix(ColSideColors))
      #stop("'ColSideColors' must be a matrix")
      if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      #lhei <- c(lhei[1], 0.2, lhei[2])
      lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
    }
    
    if (!missing(RowSideColors)) {
      #if (!is.matrix(RowSideColors))
      #stop("'RowSideColors' must be a matrix")
      if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
      #lwid <- c(lwid[1], 0.2, lwid[2])
      lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)){
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[,rowInd, drop=F])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      if (length(rownames(RowSideColors)) > 0) {
        axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  if (!missing(ColSideColors)) {
    
    if (!is.matrix(ColSideColors)){
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      par(mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[colInd, , drop=F]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, KeyValueName, line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}

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

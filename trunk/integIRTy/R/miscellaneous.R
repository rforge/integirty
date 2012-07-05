## BI computation from OOMPA Suite
bimodalIndex <- function(x){
	mc <- try(Mclust(x, G = 2, modelNames = "E"), silent=TRUE)
    if(class(mc)!='try-error'){
		sigma <- sqrt(mc$parameters$variance$sigmasq)
		delta <- abs(diff(mc$parameters$mean))/sigma
		pi <- mc$parameters$pro[1]
		bi <- delta * sqrt(pi * (1 - pi))
		res <- c(mc$parameters$mean, sigma, delta, pi, bi)
	} else {
		res <- rep(NA, 6)
	}	
    names(res) <- c("mu1", "mu2", "sigma", "delta", "pi", "BI")    
	res
}
Zscore2Binary <- function(Zscore, tau1=-2, tau2=2)
{
	binaryRes <- Zscore
	binaryRes[which(Zscore<tau2 & Zscore>tau1, arr.ind=T)] <- 0 
	binaryRes[which(Zscore>=tau2 | Zscore<=tau1, arr.ind=T)] <- 1
	binaryRes
}

# expr: tumor expression matrix
# exprCtr: normal control expression matrix, same number of rows as tumor; may contain 
#		   missing values and will propogate
# refUseMean: logical indicating whether to use mean or median of normal samples as reference.
#			default use median which is more robust
# BIthr: threshold of bimodality index to flag bimodal genes
# parallel: logical indicating whether to use parallel backend for BI computation; 
#			since the computation is fast, when gene number is small, parallel is actually slower.
dichotomizeExpr <- function(expr, exprCtr, refUseMean=FALSE, BIthr=NULL, tau1=-2.5, tau2=2.5, parallel=FALSE){
	if(is.null(BIthr)) BIthr <- ifelse(ncol(expr)>50, 1.1, 2) # guess BI threshold
	## BI computation for tumor samples
	compBI <- function(ind) {bimodalIndex(expr[ind, ])}
	if(parallel==TRUE){
		# call a function defined elsewhere; might be changed to: .packages='intIRT'
		BIinfo <- foreach(i=seq_len(nrow(expr)), .combine='rbind', .export='bimodalIndex') %dopar% compBI(i)
	} else {
		BIinfo <- foreach(i=seq_len(nrow(expr)), .combine='rbind') %do% compBI(i)
	}	
	
	## computation of Zmat
	if(refUseMean==TRUE) {
		normalRef <- apply(exprCtr, 1, mean, na.rm=T)
	} else {
		normalRef <- apply(exprCtr, 1, median, na.rm=T)
	}
	sd_normal <- apply(exprCtr, 1, sd, na.rm=T)
	sd_tumor <- apply(expr, 1, sd, na.rm=T)
	sd_tumor_2component <- BIinfo[, 'sigma'] # 2-component SD
	BIs <- BIinfo[, 'BI']
	indBimodal <- which(BIs>=BIthr) # index for bimodal genes; works even no bimodal genes are selected
	sd_tumorUsed <- sd_tumor # SD of tumor to be used on the computation
	sd_tumorUsed[indBimodal] <- sd_tumor_2component[indBimodal]
	# plot(sd_tumor, sd_tumorUsed); abline(0, 1, col=4, lty=2)
	# plot(sd_normal, sd_tumorUsed); abline(0, 1, col=4, lty=2)
	sd_adjusted <- sqrt(sd_tumorUsed^2+sd_normal^2/ncol(exprCtr)) # adjusted SD as in the denominator
	# plot(sd_tumorUsed, sd_adjusted); abline(0, 1, col=4, lty=2)
	matTumor_centered <- expr - normalRef
	Zmat <- matTumor_centered
	Zmat <- sweep(matTumor_centered, 1, sd_adjusted, '/')
	## binary mat
	BinaryMat <- Zscore2Binary(Zmat, tau1=tau1, tau2=tau2)	
	BinaryMat	
}
# no parallel is needed since all comutation are vectorized
dichotomizeMethy <- function(methy, methyCtr, refUseMean=FALSE){
	if(refUseMean==TRUE) {
		normalRef <- apply(methyCtr, 1, mean, na.rm=T)
	} else {
		normalRef <- apply(methyCtr, 1, median, na.rm=T)
	}
	myCut <- function(t, labels=c('low', 'mid', 'high'))
	{	
		## there are several instances of beta being 0, hence, the first break should be less than 0, otherwise, NA returned
		cut(t, breaks=c(-0.01, 0.25, 0.75, 1), labels=labels)
	}
	
	# mat of 3 status formed by cut
	TriMatMethyByCut <- t(apply(methy, 1, myCut))
	ref <- as.character(myCut(normalRef))
	res <- matrix(as.numeric(sweep(TriMatMethyByCut, 1, ref, '!=')), byrow=F, nrow=nrow(methy))
	rownames(res) <- rownames(methy)
	colnames(res) <- colnames(methy)
	res
}
dichotomizeCN <- function(CN, CNctr=NULL, tau1=-0.3, tau2=0.3){
	if(!is.null(CNctr) & ncol(CN)==ncol(CNctr)){
		mat <- CN-CNctr # when paired samples are specified, correct for germline CN change
	} else {
		mat <- CN # when no normal control available or not paired, no correction is performed
	}
	Zscore2Binary(mat, tau1=tau1, tau2=tau2)
}
# a wrapper of the 3 dichtomizing methods
dichotomize <- function(mat, matCtr, assayType=c('Expr', 'Methy', 'CN'), ...){
	assayType <- try(match.arg(assayType, c('Expr', 'Methy', 'CN'), several.ok=FALSE), silent=TRUE)
	# stop if assays not recognizable.
	if(class(assayType)=='try-error') stop('Assays other than expression, methylation or copy number is specified. Please
		dichotomize them mannually!\n')
	res <- switch(assayType, 'Expr'=dichotomizeExpr(mat, matCtr, ...),
					  'Methy'=dichotomizeMethy(mat, matCtr, ...),
					  'CN'=dichotomizeCN(mat, matCtr, ...))	
	res				  
}

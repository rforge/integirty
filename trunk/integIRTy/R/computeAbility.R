computeAbility <-
function(respMat, dscrmn=dscrmn, dffclt=dffclt, c=rep(0, length(dffclt)), parallel=FALSE)
{
	d=rep(1, length(dffclt))
	method='BM'
	# item parameter input matrix
	itemParMat <- cbind(dscrmn, dffclt, c, d)
	# extract unique pattern:
	PatternStr <- apply(respMat, 1, paste, collapse='')
	uniquePatternStr <- unique(PatternStr)
	unique2allInd <- match(uniquePatternStr, PatternStr)
	all2uniqueInd <- match(PatternStr, uniquePatternStr)
	uniquePatternMat <- matrix(as.numeric(unlist(strsplit(uniquePatternStr, ''))), ncol=ncol(respMat), byrow=TRUE)

	thetaEstFunc <- function(i) thetaEst(itemParMat, uniquePatternMat[i, ], method=method, range=c(-8,8))
	if(parallel==TRUE){
		abilityEst <- foreach(i=1:nrow(uniquePatternMat), .combine='c', .packages='integIRTy') %dopar% thetaEstFunc(i)
	} else {
		abilityEst <- foreach(i=1:nrow(uniquePatternMat), .combine='c') %do% thetaEstFunc(i)
	}
	names(abilityEst) <- uniquePatternStr
	
	allScore <- abilityEst[all2uniqueInd]
	names(allScore) <- rownames(respMat) # respMat contains no NA, hence rowname is corresponding gene name
	#list(uniqueScore=abilityEst, unique2allInd=unique2allInd, all2uniqueInd=all2uniqueInd, compTime=compTime, allScore=allScore)
	allScore
}


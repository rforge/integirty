calculatePermutedScoreByGeneSampling <-
function(originalMat, dscrmn=dscrmn, dffclt=dffclt, c=rep(0, length(dffclt)), fold=1, parallel=FALSE){
	nGeneSampled <- floor(nrow(originalMat)*fold)			
	permutedMat <- foreach(i=1:ncol(originalMat), .combine='cbind') %do% {
		rbinom(nGeneSampled, 1, mean(originalMat[, i]))
	}
	computeAbility(permutedMat, dscrmn=dscrmn, dffclt=dffclt, c=c, parallel=parallel)
	
}


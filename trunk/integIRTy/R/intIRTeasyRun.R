intIRTeasyRun <- 
	function(platforms, model=3, guessing=FALSE, addPermutedScore=FALSE, fold=1, echo=TRUE, parallel=FALSE){
	nPlat <- length(platforms)
	########################################
	# fit on individual platform
	########################################
	if(echo==TRUE){
			cat('Performing ltm fit for each platform', '\n')	
	}	
	#fits <- list('vector')
	fitFunc <- function(i) fitOnSinglePlat(data=platforms[[i]], model=model, guessing=guessing)	
	## fits: a list of list
	if(parallel==TRUE){
		fits <- foreach(plat = 1:nPlat, .packages='integIRTy') %dopar% fitFunc(plat)
	} else {
		fits <- foreach(plat = 1:nPlat) %do% fitFunc(plat)
	}	
	########################################
	# calculate latent trait estimates from each platform
	########################################
	if(echo==TRUE){
			cat('Performing latent trait estimation for each platform', '\n')	
	}
	# assemble the item parameters and input matrix as a list: each individual platform followed by integrated data
	a_list <- foreach(i=1:nPlat) %do% coef(fits[[i]]$fit)[, 'Dscrmn']; a_list[[length(a_list)+1]] <- c(unlist(a_list)) # append integrated parameter
	b_list <- foreach(i=1:nPlat) %do% coef(fits[[i]]$fit)[, 'Dffclt']; b_list[[length(b_list)+1]] <- c(unlist(b_list)) # append integrated parameter
	tempf <- function(i) {
		if(guessing==TRUE) {
			coef(fits[[i]]$fit)[, 'Gussng']
		} else {
			rep(0, ncol(platforms[[i]]))
		}
	}	
	c_list <- foreach(i=1:nPlat) %do% tempf(i); c_list[[length(c_list)+1]] <- c(unlist(c_list)) # append integrated parameter
	platform_list <- foreach(i=1:nPlat) %do% platforms[[i]]; 
	platform_list[[length(platform_list)+1]] <- foreach(i=1:nPlat, .combine='cbind') %do% platforms[[i]]
	compScoreFunc <- function(i){
		computeAbility(platform_list[[i]], 
						dscrmn=a_list[[i]],
						dffclt=b_list[[i]], 
						c=c_list[[i]])
	}
	# this is a matrix each column corresponding to a platform in the original order
	if(parallel==TRUE){
		estimatedScoreMat <- foreach(plat = 1:length(platform_list), .combine='cbind', .packages='integIRTy') %dopar% compScoreFunc(plat)
	} else {
		estimatedScoreMat <- foreach(plat = 1:length(platform_list), .combine='cbind') %do% compScoreFunc(plat)
	}
	
	########################################
	# if necessary, add permuted score for each platform as well as integrated data
	########################################
	if(addPermutedScore==TRUE){
		if(echo==TRUE){
			cat('Calculating permuted latent trait for each platform', '\n')	
		}	
		compPermutedScoreFunc <- function(i) {
			calculatePermutedScoreByGeneSampling(platform_list[[i]], 
						dscrmn=a_list[[i]],
						dffclt=b_list[[i]], 
						c=c_list[[i]], fold=fold)
		}
		if(parallel==TRUE){
			permutedScoreMat <- foreach(plat = 1:length(platform_list), .combine='cbind', .packages='integIRTy') %dopar% compPermutedScoreFunc(plat)
		} else {
			permutedScoreMat <- foreach(plat = 1:length(platform_list), .combine='cbind') %do% compPermutedScoreFunc(plat)
		}
		
	}
	# returns integrated latent trait and ltm fit on individual platform					
	if(addPermutedScore==TRUE){
		return(list(fits=fits, estimatedScoreMat=estimatedScoreMat, permutedScoreMat=permutedScoreMat))
	} else {
		return(list(fits=fits, estimatedScoreMat=estimatedScoreMat))
	}
}

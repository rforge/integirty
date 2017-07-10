intIRTeasyRunFromRaw <- function(platforms, platformsCtr, assayType=c('Expr', 'Methy', 'CN'),
                                 model=3, guessing=FALSE, permutationMethod=NULL, fold=1,
                                 nPerm=200, echo=TRUE, parallel=FALSE, ...) {
  assayType <- try(match.arg(assayType, c('Expr', 'Methy', 'CN'), several.ok=TRUE), silent=TRUE)
  if (is.null(permutationMethod)) permutationMethod <- "none"
	# stop if assays not recognizable.
  if(class(assayType)=='try-error') stop('Assays other than expression, methylation or copy number is specified. Please
		dichtomize them mannually and pass them to intIRTeasyRun!\n')
  nPlat <- length(platforms)
  if(nPlat!=length(platformsCtr)) stop('Number of assays in tumor samples should equal to the number in normal controls!\n')
  if(nPlat!=length(assayType)) stop('Number of assays in tumor samples should equal to the number assay types specified!\n')
	
	# stop when requiring sample label permutation but miss normal control samples
  if(sum(is.na(platformsCtr))>0 & permutationMethod=='sample label permutation')
    stop('Normal control samples for all assays are needed for sample label permutation!\n')
	# stop if permutation approach is mis-specified
  if(permutationMethod != 'none'){
    permutationMethod <- try(match.arg(permutationMethod, c('gene sampling', 'sample label permutation'),
                                       several.ok=FALSE), silent=TRUE)
    if(class(permutationMethod)=='try-error') stop('Permutation method can only be: gene sampling, sample label permutation or NULL!\n')
  }	
	########################################
	# dichtomize and easyRun on the dichotomized datasets
	########################################
  dichotomizeFunc <- function(i) {
    if(assayType[i]=='Expr') {
      res <- dichotomize(platforms[[i]], platformsCtr[[i]], assayType=assayType[i],
                         parallel=parallel) # parallel only happens on expression data
    } else {
      res <- dichotomize(platforms[[i]], platformsCtr[[i]], assayType=assayType[i])
    }
    res
  }
  if(echo==TRUE) cat('Performing data dichtomization for each platform', '\n')
	# no parallel on dichtomization part except Expr; 
  i <- NULL # to quite R CMD check
  platformBinary <- foreach(i=1:nPlat, .packages='integIRTy') %do% dichotomizeFunc(i)

  if(permutationMethod == 'none') {
    addPermutedScore <- FALSE; # no permutation when addPermutedScore=NULL
  } else if(permutationMethod=='gene sampling') {
    addPermutedScore <- TRUE; # gene sampling specified
  } else {
    addPermutedScore <- FALSE; # sample permutation specified: add further computation
  }
  if(echo==TRUE) cat('Performing easyRun on dichtomized data', '\n')
	
  easyRunResBinary <-  intIRTeasyRun(platformBinary, model=model, guessing=guessing, 
                                     addPermutedScore=addPermutedScore, fold=fold, echo=echo, parallel=parallel)
	
	########################################
	# add sample label permutation if necessary
	########################################	
	## item parameters
  a_list <- easyRunResBinary[['dscrmnList']]			
  b_list <- easyRunResBinary[['dffcltList']]			
  c_list <- easyRunResBinary[['gussngList']]		
	
  sampleLabelPermFunc <- function(i) {
		# binary matrices after sample label permutation
    plat <- NULL # to quite R CMD check
    platformPermuted <- foreach(plat=1:nPlat, .packages='integIRTy') %do% {
      tempMat <- cbind(platforms[[plat]], platformsCtr[[plat]])
      tempMat <- tempMat[, sample(1:ncol(tempMat))] # permute data
      dichotomize(tempMat[, 1:ncol(platforms[[plat]])], tempMat[, -(1:ncol(platforms[[plat]]))], assayType=assayType[plat])
    }
    platformPermuted_list <- platformPermuted; 
    platformPermuted_list[[nPlat+1]] <- foreach(i=1:nPlat, .combine='cbind') %do% platformPermuted[[plat]]
    compScoreFunc <- function(j){
      computeAbility(platformPermuted_list[[j]], 
                     dscrmn=a_list[[j]],
                     dffclt=b_list[[j]], 
                     c=c_list[[j]])
    }
		# permuted score mat at ith iteration: nPlat+1 columns with the last column for integrated data
		# no parallel within iteration
    if(i %% 10==0) cat('.')
    if(i %% 100==0) cat('\n')
    plat <- NULL # to quiet R CMD check
    tempPermutedScoreMat <- foreach(plat=1:length(platformPermuted_list), .combine='cbind') %do% compScoreFunc(plat)
    tempPermutedScoreMat
  }
  combMatToArray <- function(x, y) abind(x, y, along=3)
	# permutedScoreMatWithLabelPerm: array with the 3rd dimension to be different iterations
  i <- NULL # to quiet R CMD check
  if(permutationMethod=='sample label permutation'){
    if(echo==TRUE) cat('Performing sample label permutation', '\n')
    if(parallel==TRUE){
      permutedScoreMatWithLabelPerm <- foreach(i = 1:nPerm, .combine='combMatToArray') %dopar% sampleLabelPermFunc(i)
    } else {
      permutedScoreMatWithLabelPerm <- foreach(i = 1:nPerm, .combine='combMatToArray') %do% sampleLabelPermFunc(i)
    }
		# append the permuted result
    easyRunResBinary[['permutedScoreMatWithLabelPerm']] <- permutedScoreMatWithLabelPerm
  }	
  easyRunResBinary
}


fitOnSinglePlat <-
function(data, model=2, guessing=FALSE, 
							sampleIndices=1:ncol(data), 
							geneIndices=1:nrow(data), ...)
{
	## tpm fits a guessing parameter in 2 modes: one is guessing in unconstrained rasch,
	## the other is guessing in 2PL.
	## if guessing=T, then max.guessing=1 is used to fit the guessing parameter for every item
	#K <- length(sampleIndices) # item size
	if(is.null(sampleIndices)) sampleIndices=1:ncol(data)
	if(is.null(geneIndices)) sampleIndices=1:nrow(data)
	
	temp <- data[geneIndices, sampleIndices]
	if(model==1)
	{
		tempFit <- try(rasch(temp, constraint = cbind(ncol(temp) + 1, 1), ...), silent=TRUE)
	}
	if(model==2)
	{# force all dscrmn the same. dscmn is contained as a pure vector in the coef.
		if(guessing==FALSE) # no guessing
		{
			tempFit <- try(rasch(temp, ...), silent=TRUE)
		} else # with guessing
		{
			tempFit <- try(tpm(temp, type = "rasch", ...) , silent=TRUE)
		}
	}
	if(model==3)
	{
		if(guessing==FALSE) # no guessing
		{
			tempFit <- try(ltm(temp ~ z1, ...), silent=TRUE)
		} else # with guessing
		{
			tempFit <- try(tpm(temp, ...) , silent=TRUE)
		} 	
	}
	
	return(list(fit=tempFit, model=model, guessing=guessing, sampleIndices=sampleIndices, geneIndices=geneIndices))
}


simulateBinaryResponseMat <-
function(a=a, b=b, theta=theta){
	if(length(a)!=length(b)) stop('length of difficulty and discrimination should be equal!\n')
	nItem <- length(a)
	nExaminee <- length(theta)
	# P_ij for gene i in sample j
	P <- matrix(NA, nrow=length(theta), ncol=nItem)
	for(i in 1:nItem){
		P[, i] <- exp(a[i]*(theta-b[i]))/(1+exp(a[i]*(theta-b[i])))
	}
	# binary matrix
	X <- matrix(NA, nrow=nExaminee, ncol=nItem)
	for(i in 1:nItem){
		X[, i] <- rbinom(nExaminee, size=1, prob=P[, i])
	}
	if(!is.null(names(a))) colnames(X) <- names(a)
	if(!is.null(names(theta))) rownames(X) <- names(theta)
	X
}


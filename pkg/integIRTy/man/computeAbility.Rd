\name{computeAbility}
\alias{computeAbility}
\title{
Calculate latent traits for a given response matrix and item parameters using MLE
}
\description{
This function calculates the MLE of latent traits for a given response matrix with rows being examinees and columns being
items for given item parameters.
}
\usage{
computeAbility(respMat, dscrmn = dscrmn, dffclt = dffclt,
		c = rep(0, length(dffclt)), parallel=FALSE)
}

\arguments{
  \item{respMat}{
The response matrix of 0 and 1 with rows being examinees and columns being items.
}
  \item{dscrmn}{
A vector of item discrimination parameter.
}
  \item{dffclt}{
A vector of item difficulty parameter.
}
  \item{c}{
A vector of guessing parameter. Default is set to all 0 indicating no guessing allowed.
}
  \item{parallel}{
Logical indicating whether to use parallel computing with foreach package as backend.
}
}
\details{
This function is a wrapper of the thetaEst() function from catR package (Magis, 2012). 
}
\value{
A vector of latent trait estimates for each examinee.
}
\references{
David Magis, Gilles Raiche (2012). Random Generation of Response Patterns under Computerized Adaptive Testing with the 
R Package catR. Journal of Statistical Software, 48(8), 1-31. 
}
\author{
Pan Tong (nickytong@gmail.com), Kevin R Coombes (krc@silicovore.com)
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fitOnSinglePlat}, \link{intIRTeasyRun}, \link{calculatePermutedScoreByGeneSampling}}
}
\examples{
# number of items and number of genes
nSample <- 10
nGene <- 2000
set.seed(1000)
a <- rgamma(nSample, shape=1, scale=1)
b <- rgamma(nSample, shape=1, scale=1)
# true latent traits
theta <- rnorm(nGene, mean=0)

# probability of correct response (P_ij) for gene i in sample j
P <- matrix(NA, nrow=nGene, ncol=nSample)
for(i in 1:nSample){
	P[, i] <- exp(a[i]*(theta-b[i]))/(1+exp(a[i]*(theta-b[i])))
}
# binary matrix
X <- matrix(NA, nrow=nGene, ncol=nSample)
for(i in 1:nSample){
	X[, i] <- rbinom(nGene, size=1, prob=P[, i])
}
# IRT fitting
fit2PL <- fitOnSinglePlat(X, model=3)
dffclt <- coef(fit2PL$fit)[, 'Dffclt']
dscrmn <- coef(fit2PL$fit)[, 'Dscrmn']
# estimated latent trait
score <- computeAbility(X, dffclt=dffclt, dscrmn=dscrmn)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{regression}

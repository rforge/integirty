\name{calculatePermutedScoreByGeneSampling}
\alias{calculatePermutedScoreByGeneSampling}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the permuted latent trait by gene sampling
}
\description{
Given the original binary matrix and item parameters, calculate the permuted latent trait by
gene sampling. Basically this function permutes within columns and recompute the latent trait
using pre-specified item parameters and the permuted binary matrix. 
}
\usage{
calculatePermutedScoreByGeneSampling(originalMat, dscrmn = dscrmn,
      dffclt = dffclt, c = rep(0, length(dffclt)), fold = 1, parallel=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{originalMat}{
The original response matrix
}
  \item{dscrmn}{
The estimated item discrimination parameter
}
  \item{dffclt}{
The estimated item difficulty parameter
}
  \item{c}{
The estimated item guessing parameter if available
}
  \item{fold}{
The fold relative to the number of genes present should gene sampling achieve. Default is 1, meaning equal number
of genes are sampled. Increasing fold would increase the precesion in estimating empirical P value
}
  \item{parallel}{
Logical indicating whether to use parallel computing with foreach package as backend.
}
}
\details{
Both gene sampling and sample label permutation can be used to infer the null distribution of altent
traits. For sample label permutation, one can simply first construct the binary matrix
after permuting the sample labels and feed it to computeAbility() function together with item parameters
}
\value{
A vector of null latent traits by gene sampling
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
\code{\link{fitOnSinglePlat}, \link{intIRTeasyRun}, \link{computeAbility}}
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
# estimated null latent trait by gene sampling
scoreNull <- calculatePermutedScoreByGeneSampling(X, dffclt=dffclt,
	  dscrmn=dscrmn, fold=1)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{nonparametric}

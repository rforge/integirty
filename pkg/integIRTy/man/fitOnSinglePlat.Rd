\name{fitOnSinglePlat}
\alias{fitOnSinglePlat}
\title{
Fit IRT model on a single platform
}
\description{
This function fits the Item Response Model for one platform. It assumes the user has already dichotomized the data. 
}
\usage{
fitOnSinglePlat(data, model = 2, guessing = FALSE,
    sampleIndices = 1:ncol(data), geneIndices = 1:nrow(data), ...)
}
\arguments{
  \item{data}{
    A matrix of 0's and 1's with rows being genes (treated as examinees) and columns being samples (treated as items).
  }
  \item{model}{
    IRT model. 1-Rasch model where all item discrination are set to 1;
    2-all item discrimation are set to be equal but not necessarily as 1; 
    3-the 2PL model where no constraint is put on the item difficulty and discrimination parameter.  
}
  \item{guessing}{
A logical variable indicating whether to include guessing parameter in the model.
}
  \item{sampleIndices}{
Indices of the samples to be feeded into the model. Default is set to use all samples.
}
  \item{geneIndices}{
Indices of the genes to be feeded into the model. Default is to use all genes.
}
  \item{\dots}{
Additional options available in ltm package. Currently not used in intIRT package.
}
}

\value{
A list giving the estimated IRT model and related information
\item{fit}{An object returned by calling ltm package. Item parameters and other auxillary inforamtion (i.e. loglikelihood, convergence, Hessian) can be accessed
	from this object. For more details, please refer to ltm package}
\item{model}{The model type}
\item{guessing}{The guessing parameter}
\item{sampleIndices}{The sample indices used in the model}
\item{geneIndices}{The gene indices used in the model}
}
\references{
Rizopoulos, D. (2006) ltm: An R package for latent variable modelling and item response theory analyses. 
Journal of Statistical Software, 17(5), 1-25. URL http://www.jstatsoft.org/v17/i05/
}
\author{
Pan Tong (nickytong@gmail.com), Kevin R Coombes (krc@silicovore.com)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{computeAbility}, \link{intIRTeasyRun}, \link{calculatePermutedScoreByGeneSampling}}
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

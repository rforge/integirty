\name{simulateBinaryResponseMat}
\alias{simulateBinaryResponseMat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate binary response matrix according to 2-parameter Item Characteristic Function for
given latent traits and item parameters. 
}
\description{
This function generates binary response matrix according to the Item Characteristic Function for specified
item parameter and latent traits. It can be used for simulation purposes.
}
\usage{
simulateBinaryResponseMat(a = a, b = b, theta = theta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
A vector of item discrimination parameter 
}
  \item{b}{
A vector of item difficulty parameter 
}
  \item{theta}{
A vector of true latent traits
}
}
\details{
This function is not necessary for the integration purpose. It serves as a utility function to help the
user conduct simulation.
}
\value{
A matrix of 0's and 1's where rows are genes (examinees) and columns are samples (items).
}


\author{
Pan Tong (nickytong@gmail.com), Kevin R Coombes (kcoombes@mdanderson.org)
}



%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{computeAbility}, \link{fitOnSinglePlat}, \link{intIRTeasyRun}}
}
\examples{
# number of samples and genes to simulate
nSample <- 50
nGene <- 1000
# mean and variance of item parameters
meanDffclt_Expr <- 3; varDffclt_Expr <- 0.2
meanDscrmn_Expr <- 1.5; varDscrmn_Expr <- 0.1
# generate item parameters from gamma distribution
set.seed(1000)
Dffclt_Expr <-  rgamma(nSample, shape=meanDffclt_Expr^2/varDffclt_Expr,
                      scale=varDffclt_Expr/meanDffclt_Expr)
Dscrmn_Expr <-  rgamma(nSample, shape=meanDscrmn_Expr^2/varDscrmn_Expr,
                      scale=varDscrmn_Expr/meanDscrmn_Expr)
# generate latent trait
theta <- rnorm(nGene)
# the binary response matrix
binary_Expr <- simulateBinaryResponseMat(a=Dscrmn_Expr, b=Dffclt_Expr, theta=theta)
dim(binary_Expr)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

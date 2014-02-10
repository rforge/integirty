\name{intIRTeasyRun}
\alias{intIRTeasyRun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The easyrun function for integrating multiple modalities of high throughput assays using binary input matrix. 
}
\description{
It fits IRT models on each
of the specified platform and calculate integrated latent trait. If required, permuted latent trait by gene sampling
will also be calculated. An option for parallel computing is implemented to speed up the computation. 
}
\usage{
intIRTeasyRun(platforms, model = 3, guessing = FALSE,
    addPermutedScore = FALSE, fold = 1, echo = TRUE, parallel = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{platforms}{
A list of response matrices representing different platforms. It assumes that the number of rows (genes ) must be equal whiel the number of columns (samples)
can be different.
}
  \item{model}{
The model type as described in fitOnSinglePlat.
}
  \item{guessing}{
A logical variable indicating whether to include guessing parameter in the model.
}
  \item{addPermutedScore}{
A logical variable indicating whether to also calculate permuted latent trait by gene sampling.
}
  \item{fold}{
The fold of sampling to calculate permuted score as used in calculatePermutedScoreByGeneSampling(). Only relevant when addPermutedScore=TRUE is used.
}
  \item{echo}{
A logical variable indicating whether to print out the progress information.
}
  \item{parallel}{
Logical indicating whether to use parallel computing with foreach package as backend.
}
}
\details{
Parallel computing uses foreach and related packages for backend. The parallelism assumes computation on each 
platform individually takes similar time; the latent trait computation of the integrated data is assumed to be
 comparable to computation on individual platform. By default, all parallel options are set to be FALSE. 
Parallelism happens on the individual assay and combined data level;
No parallelism happens on genes since it would only slow donw the computation due to data transfering!
}
\value{
A list with following elements:
\item{fits }{Model fits for each platform as returned by fitOnSinglePlat function}
\item{estimatedScoreMat}{A matrix of estimated latent traits. The first several columns correspond to the
	individual assays; the last column represents the integrated latent trait with all data.}
\item{permutedScoreMat }{{A matrix of latent trait estimates after permuting the binary matrix within columns. This is only available
	if addPermutedScore is set to TRUE. The first several columns correspond to the
	individual assays; the last column represents the integrated data.}
}
\item{dscrmnList }{A list of discrimination parameters. Each element contains all of the discrimination parameters as a vector for each
assay. The last element contains the discrimination parameters for the integrated data which is formed by combining discrimination
parameters from each assay sequentially.}
\item{dffcltList }{ Same format as dscrmnList except it contains difficulty parameter.
}
\item{gussngList }{ Same format as dscrmnList except it contains guessing parameter. Be default, this is just all 0's.
}
}
\references{
Tong, P. and Coombes, K. R. (2012). integIRTy: a method to identify altered genes in cancer accounting for
multiple mechanisms of regulation using item response theory. Bioinformatics, submitted.
}
\author{
Pan Tong (nickytong@gmail.com), Kevin R Coombes (krc@silicovore.com)
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{intIRTeasyRunFromRaw}, \link{fitOnSinglePlat}, \link{calculatePermutedScoreByGeneSampling}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

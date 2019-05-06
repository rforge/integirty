\name{intIRTeasyRunFromRaw}
\alias{intIRTeasyRunFromRaw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The easyrun function for integrating multiple modalities of high throughput assays using raw data. 
}
\description{
This function performs data dichotomization, IRT fitting on individual assay, latent trait estimation for integrated data and significance
assessment of latent trait by permutation. An option for parallel computing is implemented to speed up the computation. 
}
\usage{
intIRTeasyRunFromRaw(platforms, platformsCtr, 
	assayType = c("Expr", "Methy", "CN"), 
	model = 3, guessing = FALSE, permutationMethod = NULL, 
	fold = 1, nPerm = 200, echo = TRUE, 
	parallel = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{platforms}{
A list of matrices of the raw data for tumor samples. The matrices should have equal row number corresponding to the same set of genes. The columns representing the tumor
samples can differ.
}
  \item{platformsCtr}{
A list of matrices of the raw data for normal control samples. The matrices should have equal row number corresponding to the same set of genes. The column number of each
matrix can differ. When normal control is not available, i.e. in CN data, use NA instead.
}
  \item{assayType}{
A vector of assay types. Candidates can only be a combination of "Expr", "Methy", "CN" in the order of the assays specified in the input platforms. 
When assays other then these three types, we recommend the user to dichotomize the data first and use the intIRTeasyRun function. 
}
  \item{model}{
The model type as described in fitOnSinglePlat. 1: Rasch model where all item discrination are set to 1; 2: all item discrimation are set to be equal but not necessarily as 1;
3: the 2PL model where no constraint is put on the item difficulty and discrimination parameter.  
}
  \item{guessing}{
A logical variable indicating whether to include guessing parameter in the model.
}
  \item{permutationMethod}{
What permutation method to use. It can only be 'gene sampling', 'sample label permutation' or NULL. if NULL,
no permutation is performed
}
  \item{fold}{
The fold of sampling to calculate permuted score as used in calculatePermutedScoreByGeneSampling(). 
Only relevant when permutationMethod=gene sampling is used.
}
  \item{nPerm}{
Number of permutations for sample label permutation. It is effective only when permutationMethod='sample label permutation'.
}
  \item{echo}{
A logical variable indicating whether to print out the progress information.
}
  \item{parallel}{
Logical indicating whether to use parallel computing with foreach package as backend.
}
  \item{\dots}{
Additional parameters for dichotomizing binary data.
}
}
\value{
A list quite similar to the results returned by intIRTeasyRun. The following elements are included:
\item{fits }{Model fits for each platform as returned by fitOnSinglePlat function}
\item{estimatedScoreMat}{A matrix of estimated latent traits. The first several columns correspond to the
	individual assays; the last column represents the integrated latent trait with all data.}
\item{permutedScoreMat }{{A matrix of latent trait estimates after permuting the binary matrix within columns. This is only available
	if permutationMethod='gene sampling' is used. The first several columns correspond to the
	individual assays; the last column represents the integrated data.}
}
\item{dscrmnList }{A list of discrimination parameters. Each element contains all of the discrimination parameters as a vector for each
assay. The last element contains the discrimination parameters for the integrated data which is formed by combining discrimination
parameters from each assay sequentially.}
\item{dffcltList }{ Same format as dscrmnList except it contains difficulty parameter.
}
\item{gussngList }{ Same format as dscrmnList except it contains guessing parameter. Be default, this is just all 0's.
}
\item{permutedScoreMatWithLabelPerm }{ A matrix of latent trait estimates using sample label permutation. This is only available
	if permutationMethod='sample label permutation' is used. The first several columns correspond to the
	individual assays; the last column represents the integrated data.
}
}
\references{
  Tong P, Coombes KR. \emph{integIRTy: a method to identify altered
    genes in cancer accounting for multiple mechanisms of regulation
    using item response theory}. Bioinformatics, 2012 Nov 15;
  28(22):2861--9.  
}
\author{
Pan Tong (nickytong@gmail.com), Kevin R Coombes (krc@silicovore.com)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{intIRTeasyRun}, \link{fitOnSinglePlat}, \link{calculatePermutedScoreByGeneSampling}}
}
\examples{
data(OV)
# 
controlList <- list(Expr_N, Methy_N, CN_N)
tumorList <- list(Expr_T, Methy_T, CN_T)
# not run as it takes time
#runFromRaw <- intIRTeasyRunFromRaw(platforms=tumorList, 
#		platformsCtr=controlList, 
#		assayType=c("Expr", "Methy", "CN"), 
#		permutationMethod="gene sampling")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

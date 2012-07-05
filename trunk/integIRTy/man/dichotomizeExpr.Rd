\name{dichotomizeExpr}
\alias{dichotomizeExpr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dichotomize the expression data given both tumor and normal samples.
}
\description{
This function implements the  z-like metric described in the paper.
}
\usage{
dichotomizeExpr(expr, exprCtr, refUseMean = FALSE, BIthr = NULL, tau1 = -2.5, tau2 = 2.5, parallel = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expr}{
The expression matrix for tumor samples. Rows are genes and columns are samples.
}
  \item{exprCtr}{
Expression matrix of normal controls. Genes should exactly the same as the tumor sample. The sample size are not
necessarily the same as tumor sample. 
}
  \item{refUseMean}{
Logical indicating whether to use mean of normal sample as reference. Default is set to FALSE which means to use median
as it is more robust.
}
  \item{BIthr}{
 Threshold of bimodality index to flag bimodal genes. If not specified, it will be set according to the sample size of tumor
 samples. Specifically, if tumor sample size is over 100, BIthr=1.1. If sample size is between 50 and 100, BIthr=1.5. If sample
 size is below 50, BIthr=2.0.
 }
  \item{tau1}{
Lower bound of z-like metric to be coded as 0. 
}
  \item{tau2}{
Upper bound of z-like metric to be coded as 0. The z-like metric between tau1 and tau2 will be finally converted to 0 and otherwise.
}
  \item{parallel}{
Logical indicating whether to use parallel backend provided by foreach and related packages.
 }
}
\details{
The parallelism is written to speedup BI computation. If the number of genes is not large, i.e. below 4000, we recommend not to use
parallel since this will only slow down the computation. In fact, except BI computation, all other operations are written with vector
operation.  
}
\value{
A binary matrix of the same dimension of input expr. Missing values will be propogated into binary matrix.
}
\references{
Tong, P. and Coombes, K. R. (2012). integIRTy: a method to identify altered genes in cancer accounting for
multiple mechanisms of regulation using item response theory. Bioinformatics, submitted.
}
\author{
Pan Tong (nickytong@gmail.com), Kevin R Coombes (kcoombes@mdanderson.org)
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{dichotomizeCN}, \link{dichotomizeMethy}, \link{dichotomize}}
}
\examples{
data(OV)
binDat <- dichotomizeExpr(Expr_T[1:200, ], Expr_N[1:200, ])
#binDat <- dichotomizeExpr(Expr_T[1:200, ], Expr_N[1:200, ], parallel=TRUE)
binDat[15:20, 1:2]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

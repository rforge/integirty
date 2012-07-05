\name{dichotomizeMethy}
\alias{dichotomizeMethy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dichotomize the methylation data given both tumor and normal controls.
}
\description{
This function implements the procedure for dichotomizing methylation data described
in the paper.
}
\usage{
dichotomizeMethy(methy, methyCtr, refUseMean = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{methy}{
The methylation matrix for tumor samples. Each element represents the beta value which is 
bounded between 0 and 1. Rows are genes and columns are samples.
}
  \item{methyCtr}{
Methylation matrix of normal controls. Genes should exactly the same as the tumor sample. The sample size are not
necessarily the same as tumor sample. }
  \item{refUseMean}{
Logical indicating whether to use mean of normal sample as reference. Default is set to FALSE which means to use median
as it is more robust.
}
}
\value{
A binary matrix of the same dimension of input methy. 

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
\code{\link{dichotomizeCN}, \link{dichotomizeExpr}, \link{dichotomize}}
}
\examples{
data(OV)
binDat <- dichotomizeMethy(Methy_T[1:200, ], Methy_N[1:200, ])
binDat[15:20, 1:2]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

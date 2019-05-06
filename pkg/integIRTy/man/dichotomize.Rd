\name{dichotomize}
\alias{dichotomize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A wrapper that is able to dichotomize expression, methylation and CN data
}
\description{
This function provides a common interface for the user so that data dichotomization 
can be done easily.
}
\usage{
dichotomize(mat, matCtr, assayType = c("Expr", "Methy", "CN"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
A matrix, either expression, methylation or CN
}
  \item{matCtr}{
A matrix corresponding to normal controls. 
}
  \item{assayType}{
A character string specifying the assay type. It can only be any of "Expr", "Methy", "CN". For assays none of these types,
the program will quite. To run intIRT, the user can manually dichotomize the data and feed them into intIRTeasyRun function.
}
  \item{\dots}{
Additional parameters to be passed to the specific dichotomization function.
}
}
\value{
Abinary matrix of the same dimension as input mat.

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
\seealso{
\code{\link{dichotomizeCN}, \link{dichotomizeExpr}, \link{dichotomizeMethy}}
}
\examples{
data(OV)
binDat_expr <- dichotomize(Expr_T[1:20, ], Expr_N[1:20, ], assayType='Expr')
binDat_methy <- dichotomize(Methy_T[1:20, ], Methy_N[1:20, ], assayType='Methy')
}
\keyword{manip}

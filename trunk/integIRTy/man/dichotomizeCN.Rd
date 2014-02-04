\name{dichotomizeCN}
\alias{dichotomizeCN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dichotomizing copy number data based on segmented data (i.e. log2ratio).
}
\description{
A simple dichotomization procedure is implemented for CN data that only requires
two cutoffs. 
}
\usage{
dichotomizeCN(CN, CNctr = NULL, tau1 = -0.3, tau2 = 0.3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CN}{
A matrix of gene wise copy number data for tumor samples. Rows are the genes; columns are the samples.
}
  \item{CNctr}{
A matrix of copy number data for normal samples. The program first guess if this is paired data by comparing if 
tumor and normal sample has equal sample size. If TRUE, then normal samples will be subtracted element by element to
 correct for germline CN change. Otherwise, no correction is performed. The program proceeds with only tumor data.
}
  \item{tau1}{
The lower bound for log2Ratio when converting to a code as 0.
}
  \item{tau2}{
The upper bound for log2Ratio when converting to a code as 0. Log2ratio between tau1 and tau2 is converted to 0 and 1 otherwise.
}
}
\value{
Returns a binary matrix of the same dimension. Missing value would be propogated into the binary matrix.
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
\code{\link{dichotomizeExpr}, \link{dichotomizeMethy}
}
}
\examples{
data(OV)
binDat <- dichotomizeCN(CN_T[1:20, ], CN_N[1:20, ])
binDat[15:20, 1:2]
}


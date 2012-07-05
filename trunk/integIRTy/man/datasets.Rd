\name{OV}
\alias{OV}
\alias{CN_N}
\alias{CN_T}
\alias{Methy_N}
\alias{Methy_T}
\alias{Expr_N}
\alias{Expr_T}
\docType{data}
\title{Ovarian Cancer Datasets}
\description{
  Six matrices containing a subset of TCGA ovarian cancer data.
}
\usage{
data(OV)
}
\format{
  Each of the six objects (CN_N, CN_T, Methy_N, Methy_T, Expr_N, Expr_T)
  is a matrix with rows for 1000 (matched) genes and columns as samples. Gene expression, methylation
  and copy number data for 30 tumor samples and around 10 normal samples are provided.
}
\source{
  This data is a subset of the TCGA ovarian cancer datasets. The full datasets can be downloaded through the
  TCGA data portal at: http://cancergenome.nih.gov/
}
\references{
Cancer Genome Atlas Research Network (2011). Integrated genomic analyses of ovarian carcinoma. Nature 474, 609-615
}
\keyword{datasets}


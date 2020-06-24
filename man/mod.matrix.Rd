\name{mod.matrix}
\alias{mod.matrix}
\title{
Transform model matrix in a pam
}
\description{
This is an internal function of package \code{pamfe} based on the \code{\link[plm]{model.matrix.pFormula}} function from package \code{plm}. It enables the specified transformation of the model matrix in a \code{pam}, currently only a first-difference transformation is supported. Not normally called directly.
}
\usage{
mod.matrix(object, data, model = c("within", "fd"), ...)
}
\arguments{
  \item{object}{An object of class \code{\link[plm]{pFormula}}.}
  \item{data}{A \code{\link[plm]{pdata.frame}}.}
  \item{model}{Estimation method of fixed effects panel data model. Currently only first differences ("\code{fd}") are supported.}
  \item{\dots}{Further arguments.}
}

\value{
A model matrix after the specified transformation.
}

\author{
Peter Puetz \email{ppuetz@uni-goettingen.de}
}



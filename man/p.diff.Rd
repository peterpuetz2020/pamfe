\name{p.diff}
\alias{p.diff}
\title{
Transform model matrix in a pam by first differences
}
\description{
This is an internal function of package \code{pamfe}.  It enables the first-difference transformation of model matrices in a pam. Not normally called directly.
}
\usage{
p.diff(x, cond)
}
\arguments{
  \item{x}{Model matrix of class \code{\link[plm]{pdata.frame}}.}
  \item{cond}{Individual dimension of model matrix of class \code{\link[plm]{pdata.frame}}.}
}
\value{
A model matrix after building first differences.
}
\author{
Peter Puetz \email{ppuetz@uni-goettingen.de}
}


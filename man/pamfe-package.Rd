\name{pamfe-package}
\alias{pamfe-package}
\alias{pamfe}
\docType{package}
\title{
Additive Panel Data Models with Fixed Effects
}
\description{
Fits additive panel data models with fixed effects via a first-difference transformation and estimates simultaneous confidence bands.
}
\details{
As described in Puetz and Kneib (2016), the mixed model representation of penalized splines is used as it enables both smoothing parameter estimation and inferences via simultaneous confidence bands. So far, only the first-difference transformation for removing fixed effects from the corresponding panel data model is supported. Further, only gaussian responses, homoscedastic error terms without autocorrelation are supported. The estimation of additive panel data models with fixed effects (\code{\link{pam}}) is based on the \code{\link[mgcv]{gam}} function from package \code{mgcv} and the \code{\link[plm]{plm}} function from package \code{plm}. The estimation of simultaneous confidence bands (\code{\link{scb}}) uses the \code{\link[AdaptFitOS]{scbM}} function from package \code{AdaptFitOS}. For details of the underlying statistical concepts, see Puetz and Kneib (2016), Wood (2006) and Wiesenfarth et. al (2012).
}

\author{
Peter Puetz

Maintainer: Peter Puetz <ppuetz@uni-goettingen.de>
}
\references{
Puetz, P., Kneib, T. (2016). A Penalized Spline Estimator For Fixed Effects Panel Data Models.
\url{https://www.uni-goettingen.de/de/Puetz_03_2016/534166.html}


Wiesenfarth, M., Krivobokova, T., Klasen, S., Sperlich, S. (2012).  
Direct Simultaneous Inference in Additive Models and its Application to Model Undernutrition.
Journal of the American Statistical Association, 107(500): 1286-1296.

Wood, S.N. (2006). Generalized Additive Models: An Introduction with R. Chapman & Hall/CRC.
}



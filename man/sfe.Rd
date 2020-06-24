\name{sfe}
\alias{sfe}
\title{Defining smooth terms in PAM formulae}
\description{Function used in definition of smooth terms within
 \code{pam} model formulae based on the \code{\link[mgcv]{s}} function from package \code{mgcv}. The function does not evaluate a (B-spline)
 smooth - it exists purely to help set up a model using B-spline based smooths. Note that this function only supports B-spline based smooths corresponding to the \code{"ps"} 
class for the
 smoothing basis argument \code{bs} argument of the \code{\link[mgcv]{s}} function from package \code{mgcv}.}

\usage{
sfe(..., k = -1, m = NA, sp = NULL)
}
\arguments{
  \item{\dots}{Name of the covariate that this
 smooth is a function of.}
  \item{k}{The dimension of the basis used to represent the smooth term.
 See \code{\link[mgcv]{choose.k}} for further information.}
  \item{m}{A 2 item array giving the basis and penalty order for the term  separately, defaults to c(3,2), i.e., a B-spline of degree 3 and a second order (difference) penalty.}


  \item{sp}{Any supplied smoothing parameter for this term. A positive or zero element is taken as fixed 
smoothing parameter.               A negative element signals auto-initialization. Overrides values supplied in 
\code{sp} argument to \code{\link{pam}.}}
}

\value{An object similar to the \code{\link[mgcv]{smooth.construct.ps.smooth.spec}} object from package \code{mgcv} (but with fewer output values). These \code{\link[mgcv]{smooth.construct.ps.smooth.spec}} objects define smooths and are turned into 
bases and penalties by \code{\link[mgcv]{smooth.construct}}  method functions from package \code{mgcv}. 

The returned object contains the following items:


\item{term}{A text string giving the name of the covariate that 
the term is a function of.}
  
\item{bs.dim}{The dimension of the basis used to represent the smooth.}
\item{dim}{The dimension of the smoother - i.e. the number of
 covariates that it is a function of, which is always 1.}
  
\item{label}{A suitable text label for this smooth term.}


\item{sp}{array of smoothing parameters for the term (negative for 
auto-estimation). Defaults to \code{NULL}.}
}

\author{
Peter Puetz \email{ppuetz@uni-goettingen.de}
}

\examples{
# data generation: additive model with time constant indivdual fixed effects
id <- rep(1:50,each = 10)
years <- rep(1:10,50)
x1 <- runif(500)
x2 <- runif(500)
f1 <- sin(2 * pi * (x1 - 0.5)) ^ 2
f2 <- x2 * (1 - x2)
f1_s <- f1 / sd(f1)
f2_s <- f2 / sd(f2)
fe <- rep(sample(1:100,50),each = 10)
y <- fe + f1_s + f2_s + rnorm(500,sd = 0.5)
data <- as.data.frame(cbind(id,years,y,x1,x2))

# transform data set to panel data set from type "pdata.frame" from package "plm"
pdata <- pdata.frame (data, index = c("id", "years"),
                      row.names = TRUE)

# run first-difference penalized spline panel data model with generous amount of knots
mod <- pam(y ~ sfe(x1,k = 40) + sfe(x2,k = 40),data = pdata)

# particular look at smoothing term definition
sfe(x1,k = 40)
}


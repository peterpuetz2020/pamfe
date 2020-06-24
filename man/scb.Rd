\name{scb}
\alias{scb}
\title{
Calculate simultaneous confidence bands for a fitted pam object
}
\description{
Calculates simultaneous confidence bands for the mixed model representation of penalized
splines based on the \code{\link[AdaptFitOS]{scbM}} function from package \code{AdaptFitOS}. In contrast to the 
\code{\link[AdaptFitOS]{scbM}} function, confidence bands for fixed effects panel data models using the
first-difference estimator (see Puetz and Kneib, 2016) are calculated.
}
\usage{
scb(object, drv = 0, level = 0.95, pred = 500, div = 1000, pages = 0)
}
\arguments{
  \item{object}{A \code{pam} object.}
  \item{drv}{The derivative order. Defaults to 0, i.e. the estimated function itself is plotted. First and second derivatives are supported.}
  \item{level}{Level of confidence, default is 0.95.}
  \item{pred}{Number of grid points used for the plot of confidence bands, default is 500.}
  \item{div}{Precision for the integral used for calculation of the length of the curve, default is 1000.}
  \item{pages}{The number of pages over which to spread the output. If pages=1 then all terms will be plotted on one page in an automatic fashion. If pages=0 (default) all graphics settings are left as they are.}
}
\details{The resulting confidence bands have (approximate) frequentist coverage probabilities with automatic bias correction (see references).
Makes use of the volume-of-tube formula and the corresponding code from the \code{libtube} library by Catherine Loader (see package \code{\link[locfit]{locfit}}).
}
\value{
\item{crit}{A list of critical values.}
\item{seqx}{A list of grid points for the plot of confidence bands.}
\item{Stdev}{The standard deviations of estimates.}
\item{sigma2}{The variance of the residuals.}
\item{drv}{The derivative order.}
\item{fitted}{A list of fitted values.}
\item{lcb}{A list of lower bounds of confidence bands.}
\item{ucb}{A list of upper bounds of confidence bands.}
}
\references{
Krivobokova, T., Kneib, T., and Claeskens, G. (2010).
Simultaneous confidence bands for penalized spline estimators.
Journal of the American Statistical Association, 105(490):852-863.

Puetz, P., Kneib, T. (2016). A Penalized Spline Estimator For Fixed Effects Panel Data Models.
\url{https://www.uni-goettingen.de/de/Puetz_03_2016/534166.html}


Wiesenfarth, M., Krivobokova, T., Klasen, S., Sperlich, S. (2012).
Direct Simultaneous Inference in Additive Models and its Application to Model Undernutrition.
Journal of the American Statistical Association, 107(500): 1286-1296.
}
\author{
Peter Puetz \email{ppuetz@uni-goettingen.de}
}
\seealso{\code{\link{pam}}}

\examples{
library(pamfe)
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
summary(mod)

# look at 99%-confidence bands for estimated functions and its derivatives
scb_mod <- scb(mod,level = 0.99,pages = 1)
scb_der <- scb(mod,level = 0.99,pages = 1,drv = 1)
}

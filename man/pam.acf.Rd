\name{pam.acf}
\alias{pam.acf}
\title{Compute individual-wise (partial) autocorrelation function for residuals from a pam}
\description{
Computes individual-wise autocorrelation and partial autocorrelation functions for residuals from a \code{pam} and summarizes them over all individuals.
}
\usage{
pam.acf(object)
}
\arguments{
  \item{object}{A \code{pam} object.}
}
\details{
Detected error autocorrelation patterns in panel data models can be used in a subsequent reestimation of the model which accounts for these patterns, see the example below. This may be essential to obtain valid simultaneous confidence bands in a subsequent step.
}
\value{
\item{id_acf}{A matrix with one row for the autocorrelation coefficients of the residuals for each individual.}
\item{id_pacf}{A matrix with one row for the partial autocorrelation coefficients of the residuals for each individual.}
\item{summary_acf}{A matrix containing a summary for each autocorrelation coefficient, summarized over all individuals.}
\item{summary_pacf}{A matrix containing a summary for each partial autocorrelation coefficient, summarized over all individuals.}
}
\author{
Peter Puetz \email{ppuetz@uni-goettingen.de}
}
\seealso{
\code{\link{pam}}
}

\examples{
# data generation: additive model with time constant indivdual fixed effects
set.seed(2017)
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
# and no GLS-type correction which induces a negative correlation of -0.5 in adjacent
# first-differenced error
mod <- pam(y ~ sfe(x1,k = 40) + sfe(x2,k = 40), gls=FALSE, data = pdata)

# the negative correlation can be seen in the MA(1) and AR(1) coefficients:
pam_acf<-pam.acf(mod)

# in order to account for a (known or estimated) correlation structure, we could use the correlation
# classes ("corClasses") from the "nlme" package 
library(nlme)

# read ids of the first-differenced data
index_data <- as.data.frame(mod$index_diffdata)
colnames(index_data) <- "id"

# generate ma (1) process leading to the known error correlation of (almost)
# -0.5 between adjacent errors
cs2 <- corARMA(-.999,form=~1|id,p=0,q=1)
cs2 <- Initialize(cs2,data=index_data)
# gives us id-wise correlation matrix
corMatrix(cs2)

# accordingly, fitting the model with this correlation matrix gives the same results like
# running the default gls type approach which precisely accounts for the MA(1) process with 
# MA(1) coefficient -0.5 induced by first differencing
mod1 <- pam(y ~ sfe(x1) + sfe(x2), gls=FALSE, corMatrix = corMatrix(cs2), data = pdata)
mod2 <- pam(y ~ sfe(x1) + sfe(x2), gls=TRUE, data = pdata)
summary(mod1)
summary(mod2)

# of course: in a real world data analysis, we would try to model an ARMA process
# which is indicated by the estimated AR and MA coefficients

# similarly, the homoscedasticity assumption could be scrutinized, e.g. id-wise homoscedasticity
aggregate(mod$residuals,by=list(mod$index_diffdata),FUN=sd)
# ...
}
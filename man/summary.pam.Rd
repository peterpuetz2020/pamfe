\name{summary.pam}
\alias{summary.pam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary for a PAM fit}

\description{
Takes a fitted \code{pam} object produced by \code{pam()} and produces various useful
summaries from it. The code is based on the based on \code{\link[mgcv]{summary.gam}} function from package \code{mgcv}.
}
\usage{
\method{summary}{pam}(object, freq = FALSE, ...)
}
\arguments{

\item{object}{A fitted \code{pam} object as produced by \code{pam()}.}
\item{freq}{By default p-values for parametric terms are calculated using the Bayesian estimated
 covariance matrix of the parameter estimators. If this is set to \code{TRUE} then
 the frequentist covariance matrix of the parameters is used instead.}
\item{\dots}{Other arguments.}
}


\details{ Model degrees of freedom are taken as the trace of the influence (or
hat) matrix \eqn{ {\bf A}}{A} for the model fit. 
Residual degrees of freedom are taken as number of data minus model degrees of
 freedom. 
Let \eqn{ {\bf P}_i}{P_i} be the matrix 
giving the parameters of the ith smooth when applied to the data (or pseudodata in the generalized case) and let \eqn{ {\bf X}}{X} 
be the design matrix of the model. Then \eqn{ tr({\bf XP}_i )}{tr(XP_i)} is the edf for the ith term. Clearly this 
definition causes the edf's to add up properly! An alternative version of EDF is more appropriate for p-value computation, and is based on the trace of  \eqn{ 2{\bf A} - {\bf AA}}{2A - AA}.

P-values for smooth terms are usually based on a 
test statistic motivated by an extension of Nychka's (1988) analysis of the frequentist properties 
of Bayesian confidence intervals for smooths (Marra and Wood, 2012). 
These have better frequentist performance (in terms of power and distribution under the null) 
than the alternative strictly frequentist approximation. When the Bayesian intervals have good 
across the function properties then the p-values have close to the correct null distribution 
and reasonable power (but there are no optimality results for the power). Full details are in Wood (2013), 
although what is computed is actually a slight variant in which the components of the test statistic are weighted by the iterative fitting weights.


All p-values are computed without considering uncertainty in the smoothing parameter estimates. 

In simulations the p-values have best behaviour under ML smoothness selection, with REML coming second. In general the p-values behave well, but neglecting smoothing parameter uncertainty means that they may be somewhat too low when smoothing parameters are highly uncertain. High uncertainty happens in particular when smoothing parameters are poorly identified, which can occur with nested smooths or highly correlated covariates (high concurvity).

 

By default the p-values for parametric model terms are also based on Wald tests using the Bayesian 
covariance matrix for the coefficients.  
}


\value{\code{summary.pam} produces a list of summary information for a fitted \code{pam} object. 


\item{p.coeff}{Is an array of estimates of the strictly parametric model coefficients.}


\item{p.t}{Is an array of the \code{p.coeff}'s divided by their standard errors.}


\item{p.pv}{Is an array of p-values for the null hypothesis that the corresponding parameter is zero. 
Calculated with reference to the standard normal distribution.}


\item{m}{The number of smooth terms in the model.}


\item{chi.sq}{An array of test statistics for assessing the significance of
model smooth terms. See details.}


\item{s.pv}{An array of approximate p-values for the null hypotheses that each 
smooth term is zero. Be warned, these are only approximate.}


\item{se}{Array of standard error estimates for all parameter estimates.}


\item{r.sq}{The adjusted r-squared for the model. Defined as the proportion of variance explained, where original variance and 
residual variance are both estimated using unbiased estimators. This quantity can be negative if your model is worse than a one 
parameter constant model, and can be higher for the smaller of two nested models!}



\item{edf}{Array of estimated degrees of freedom for the model terms.}


\item{residual.df}{Estimated residual degrees of freedom.}


\item{n}{Number of data used for fitting process (after applying the first-difference transformation.}


\item{np}{Number of model coefficients (regression coefficients, not smoothing parameters or other parameters of likelihood).}




\item{method}{The smoothing selection criterion used.}


\item{sp.criterion}{The minimized value of the smoothness selection criterion. What is reported is the negative log marginal likelihood or negative log restricted likelihood depending on the estimation method.}


\item{scale}{Estimated (or given) scale parameter.}


\item{family}{The family (always gaussian) and link function (always identity link) used.}


\item{formula}{The original PAM formula.}




\item{pTerms.df}{The degrees of freedom associated with each parametric term
.}


\item{pTerms.chi.sq}{A Wald statistic for testing the null hypothesis that 
each parametric term is zero.}


\item{pTerms.pv}{P-values associated with the tests that each term is 
zero. The reference distribution 
is an appropriate chi-squared when the 
scale parameter is known, and is based on an F when it is not.}


\item{cov.unscaled}{The estimated covariance matrix of the parameters (or
estimators if \code{freq=TRUE}), divided 
by scale parameter.}


\item{cov.scaled}{The estimated covariance matrix of the parameters 
(estimators if \code{freq=TRUE}).}


\item{p.table}{Significance table for parameters.}


\item{s.table}{Significance table for smooths.}


\item{p.Terms}{Significance table for parametric model terms.}
}
\references{



Marra, G and S.N. Wood (2012). Coverage Properties of Confidence Intervals for Generalized Additive 
Model Components. Scandinavian Journal of Statistics, 39(1), 53-74.

 
Nychka (1988). Bayesian Confidence Intervals for Smoothing Splines. 
Journal of the American Statistical Association 83:1134-1143.


Wood, S.N. (2013). On p-values for smooth components of an extended generalized additive model. Biometrika 100:221-228.
}
\author{
Peter Puetz \email{ppuetz@uni-goettingen.de}
}
\seealso{\code{\link{pam}}}
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
summary(mod)
}
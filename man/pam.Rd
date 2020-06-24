\name{pam}
\alias{pam}
\title{Fit an additive model for panel data with fixed effects}
\description{
Fits additive panel data models with fixed effects based on the \code{\link[mgcv]{gam}} function from package \code{mgcv} and the \code{\link[plm]{plm}} function from package \code{plm}. 
Nonparametric model components are represented by penalized B-splines with smoothing parameters selected by ML or REML. For more details see \code{\link[mgcv]{gam}} from package \code{mgcv}.
}
\usage{
pam(formula, data = list(), weights = NULL, method = "REML", knots = NULL, 
optimizer = c("outer", "newton"), control = list(), sp = NULL, gls = TRUE,
 corMatrix = list(), ...)
}
\arguments{
  \item{formula}{A pam formula which is similar to the formula for a \code{\link[stats]{lm}} except that nonparametric terms via \code{\link{sfe}} can be added to the right hand side. Note that an intercept is never provided. }
  \item{data}{A data frame of class \code{\link{pdata.frame}}.}
  \item{weights}{An optional vector of weights to be used in the fitting process. Should be NULL or a numeric vector. If non-NULL, the overall magnitude of the log likelihood is not changed, i.e. the weights are normalized (\code{weights <- weights/mean(weights)}).}
  \item{method}{The smoothing parameter estimation method. \code{"REML"} for REML estimation, including of unknown scale, \code{"P-REML"} for REML estimation, but using a Pearson estimate of the scale. \code{"ML"} and \code{"P-ML"} are similar, but using maximum likelihood in place of REML. \code{"REML"} is the default.}
  \item{knots}{This is an optional list containing user specified knot values to be used for basis construction. The user simply supplies the knots to be used, which must match up with the \code{k} value supplied (note that the number of knots is not always just \code{k}).}
  \item{optimizer}{An array specifying the numerical optimization method to use to optimize the smoothing 
parameter estimation criterion (given by \code{method}). \code{"perf"} for performance iteration. \code{"outer"} 
for the more stable direct approach. \code{"outer"} can use several alternative optimizers, specified in the 
second element of \code{optimizer}: \code{"newton"} (default), \code{"bfgs"}, \code{"optim"}, \code{"nlm"} 
and \code{"nlm.fd"} (the latter is based entirely on finite differenced derivatives and is very slow).}
  \item{control}{A list of fit control parameters to replace defaults returned by \code{\link[mgcv]{gam.control}}. Values not set assume default values.}
  \item{sp}{A vector of smoothing parameters can be provided here. Smoothing parameters must be supplied in the order that the smooth terms appear in the model formula. Negative elements indicate that the parameter should be estimated, and hence a mixture of fixed and estimated parameters is possible. If smooths share smoothing parameters then \code{length(sp)} must correspond to the number of underlying smoothing parameters.}
  \item{gls}{If this argument is TRUE (the default value), then serial error correlation inherent to the first-difference transformation to remove fixed effects is accounted for via a generalized least squares approach.}
\item{corMatrix}{This is an optional list containing matrices describing the within-individual variance and correlation structure for the errors of each individual. Such matrices can easily be generated with the help of \code{\link[nlme]{corMatrix.corStruct}}. The matrices are then used via a generalized least squares approach to account for the respective error structure. For detecting specific error structures from residual checking, see \code{pam.acf}.}
  \item{\dots}{Further arguments for passing on e.g. to \code{\link[mgcv]{gam.fit}} from package \code{mgcv} which is used for the fitting process.}
}
\details{
An additive panel data models with fixed effects is a model which is capable to include individual-specific time constant effects, nonparametric effects and strictly parametric effects jointly. The fixed effects are removed by building first differences over time. The resulting dependence structure can be accounted for via a generalized least squares approach. Nonparametric effects are represented by penalized B-splines. The tradeoff between penalizing wiggliness and penalizing badness of fit is steered by associated smoothing parameters which are estimated by (restricted) maximum likelihood. For further information, see Puetz and Kneib (2016).

Note that \code{\link[mgcv]{gam}} from package \code{mgcv} is more comprehensive (e.g. it allows for generalized additive models) and offers more options to specify. The major difference is that the \code{mgcv} package is designed for cross-sectional data and panel data models with random effects. 

Details of the default underlying fitting methods are given in Wood (2011 and 2004). A concise introduction to generalized additive models and their implementation in R is given by Wood (2006). 
}

\value{An object of class \code{pam}, similar to a \code{\link[mgcv]{gam}} object from package \code{mgcv}. A \code{pam} object has has the following elements:
\item{aic}{AIC of the fitted model: bear in mind that the degrees of freedom
 used to calculate this are the effective degrees of freedom of the model, and
 the likelihood is evaluated at the maximum of the penalized likelihood in most
 cases, not at the MLE.}



\item{assign}{Array whose elements indicate which model term (listed in
\code{pterms}) each parameter relates to: applies only to non-smooth terms.}
\item{boundary}{Did parameters end up at boundary of parameter space?}
\item{coefficients}{The coefficients of the fitted model. Parametric
 coefficients are  first, followed  by coefficients for each
 spline term in turn.}
\item{control}{The \code{gam} control list used in the fit.}



\item{converged}{Indicates whether or not the iterative fitting method converged.}
\item{db.drho}{Matrix of first derivatives of model coefficients w.r.t. log smoothing parameters.}
\item{df.null}{Null degrees of freedom.}

\item{df.residual}{Effective residual degrees of freedom of the model.}
\item{edf}{Estimated degrees of freedom for each model parameter. Penalization 
means that many of these are less than 1.}
\item{edf1}{Similar, but using alternative estimate of EDF. Useful for testing.}
\item{edf2}{This edf accounts for smoothing parameter
 uncertainty. \code{edf1} is a heuristic upper bound for \code{edf2}.}


\item{family}{Family object specifying distribution (always gaussian) and link (always identity link) used.}
\item{fitted.values}{The fitted values for the model. Note that the model is fitted on data transformed by first differences.}
\item{formula}{The model formula.}



\item{gcv.ubre}{The minimized smoothing parameter selection score: negative log marginal 
likelihood or negative log restricted likelihood.}
\item{gls}{TRUE if serial error correlation inherent to the first-difference transformation to remove fixed effects was accounted for via a generalized least squares approach.}
\item{hat}{Array of elements from the leading diagonal of the `hat' (or `influence') matrix. 
Same length as response data vector.}


\item{index_data}{The individual dimension and the time dimension of the original panel data set.}
\item{index_diffdata}{The individual dimension (the ids) of the first-differenced data set.}
\item{iter}{How many iterations were required to find the smoothing parameters?
}
\item{method}{One of \code{"REML"}, \code{"P-REML"}, \code{"ML"},
 \code{"P-ML"}, depending on the fitting 
criterion used.}
\item{model}{Model frame containing all variables needed in original model fit.}

\item{n}{Number of observation used for the fittind process, i.e. after the first-difference transformation.}


\item{nsdf}{Number of parametric, non-smooth, model terms.}
\item{optimizer}{\code{optimizer} argument to \code{pam}.}
\item{outer.info}{If `outer' iteration has been used to fit the model (see
 \code{pam} argument \code{optimizer}) then this is present and contains whatever was 
returned by the optimization routine used (currently \code{\link{nlm}} or \code{\link{optim}}). }
\item{prior.weights}{Prior weights on observations.}


\item{pterms}{\code{terms} object for strictly parametric part of model.}
\item{R}{Factor R from QR decomposition of weighted model matrix, unpivoted to be in
 same column order as model matrix (so need not be upper triangular).}
\item{rank}{Apparent rank of fitted model.}
\item{reml.scale}{The scale (RE)ML scale parameter estimate.}


\item{residuals}{The residuals for the fitted model. Note that the model is fitted on data transformed by first differences.}
\item{rV}{If present, \code{rV\%*\%t(rV)*sig2} gives the estimated Bayesian covariance matrix.}
\item{scale}{When present, the scale (as \code{sig2}).}
\item{scale.estimated}{ \code{TRUE} if the scale parameter was estimated, \code{FALSE} otherwise.}
\item{sig2}{Estimated or supplied variance/scale parameter.}

\item{smooth}{List of smooth objects, containing the basis information for each term in the 
model formula in the order in which they appear.}
\item{sp}{Estimated smoothing parameters for the model. These are the underlying smoothing 
parameters, subject to optimization.}
\item{terms}{\code{terms} object of \code{model} model frame.}
\item{Vc}{Under ML or REML smoothing parameter estimation it is possible to correct the covariance
matrix \code{Vp} for smoothing parameter uncertainty. This is the corrected version.
}
\item{Ve}{Frequentist estimated covariance matrix for the parameter
 estimators. Particularly useful for testing whether terms are zero. Not so 
useful for CI's as smooths are usually biased.}


\item{Vp}{Estimated covariance matrix for the parameters. This is a Bayesian 
posterior covariance matrix that results from adopting a particular Bayesian 
model of the smoothing process.} 
\item{weights}{Final weights used in IRLS iteration.}
\item{y}{Response data used in the fitting process, i.e. after the first-difference transformation.

}
}
\references{
Puetz, P., Kneib, T. (2016). A Penalized Spline Estimator For Fixed Effects Panel Data Models.
\url{https://www.uni-goettingen.de/de/Puetz_03_2016/534166.html}


Wood, S.N. (2011). Fast stable restricted maximum likelihood 
and marginal likelihood estimation of semiparametric generalized linear 
models. Journal of the Royal Statistical Society (B) 73(1):3-36

.
Wood, S.N. (2004). Stable and efficient multiple smoothing parameter estimation for
generalized additive models. J. Amer. Statist. Ass. 99:673-686.


Wood S.N. (2006). Generalized Additive Models: An Introduction with R. Chapman
and Hall/CRC.
}
\author{
Peter Puetz \email{ppuetz@uni-goettingen.de}
}
\seealso{
\code{\link{sfe}}, \code{\link{summary.pam}}
}

\examples{
# data generation: additive model with time constant indivdual fixed effects
library(pamfe)
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

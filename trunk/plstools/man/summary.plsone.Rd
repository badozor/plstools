\name{summary.plsone}
\alias{summary.plsone}
\alias{print.summary.plsone}
\title{Summarizing Partial Least Squares regression}
\description{
This function provides additional elements to evaluate the goodness of fit
and the quality of the model (object of class inheriting from 'plsone').
}
\usage{
\method{summary}{plsone}(object, ...)
\method{print}{summary.plsone}(x,digits = max(3, getOption("digits") - 3),...)
}
\arguments{
  \item{object}{an object of class inheriting from 'plsone'.}
  \item{x}{an object of class inheriting from 'summary.plsone'.}
  \item{digits}{the number of significant digits to use when printing.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
 decrire le calcul des elements ...\cr

+ additional details on the computation of each elemnts\cr
+ At the moment this function is very slow to compute the diagnostic elements.\cr

FT2 = nf * (n * n - 1) * qf(0.95, nf, n - nf)/(n * (n - nf))\cr
DcritX and  DcritY = \cr
}
\value{
The function 'summary.plsone' retruns an list of class 'summary.plsone' containing
several elements associated with the model diagnostic:
  \item{rx2}{a matrix containing the R2 between x and th}
  \item{ry2}{a matrix containing the R2 between y and th}
  \item{T2}{a numeric vector containing Hotelling's T2 for each observation.}
  \item{Ts2}{a numeric vector containing adapted T2 for each observation (statistic used in SIMCA-P).}
  \item{FT2}{a numerical vector containing threshold to detected outlier with Hotelling's T2 (F-value).}
  \item{rxt2}{a matrix containing the correlations between x and th}
  \item{ryt2}{a matrix containing the correlations between y and th}
  \item{DModX}{a numeric vector containing distance from the \eqn{i^{th}} observation to the model in the space of the explanatory variables}
  \item{DModXN}{a numeric vector containing normaliszed values of DModX.}
  \item{DcritX}{a numerical vector containing threshold to detected the observations badly rebuilt by the model.}
  \item{DModY}{a numeric vector containing residuals of the regression of y on the 'nf' components.}
  \item{DModYN}{a numeric vector containing normaliszed values of DModY.}
  \item{DcritY}{a numerical vector containing threshold to detected the observations badly rebuilt by the model.}
  \item{coefficients}{a numeric vector containing regression coefficients of y on the components 'th'}
  \item{r.squared}{a numerical value corresponding to the R squared of linear model
  ('fraction of variance explained by the model',see 'summary.lm').}
  \item{adj.r.squared}{a numerical value corresponding to the adjusted R squared of linear model}
  \item{bh}{an numeric vector which contains coefficients associated with the explanatory variables}
  \item{n}{observation number}
  \item{p}{number of explanatory varaibles.}
  \item{nf}{an integer indicating the number of kept components (h).}
}
\references{
Tenenhaus M. (1998) La Regression PLS. Theorie et pratique. Technip, Paris.\cr}
\seealso{\code{\link{plsone}},\code{\link[stats:summary.lm]{summary.lm}}}
\examples{
require(pls)
data(yarn)
plstest1 <- plsone(density ~ NIR, nf=6, data = yarn,scale=FALSE)
summary(plstest1)
summary(plstest1$lm)
}



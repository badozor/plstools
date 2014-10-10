\name{plsone}
\alias{plsone}
\alias{plsone.fit}
\alias{cplsone.fit}
\alias{print.plsone}
\alias{coef.plsone}
\alias{coefficients.plsone}
\alias{fitted.plsone}
\alias{residuals.plsone}
\alias{formula.plsone}
\alias{delete.intercept}
\title{PLS regression based on NIPALS algorithm}
\description{
Partial least squares regresssion based on NIPALS algorithm.
}
\usage{
plsone(formula, data,subset,na.action,nf = 2,scale=TRUE,method="plsone.fit",...)
plsone.fit(y, X, nf = 2,scale=TRUE, ...)
cplsone.fit(y, X, nf = 2,scale=TRUE, ...)
\method{print}{plsone}(x, ...)
\method{coef}{plsone}(object,type="components", ...)
\method{coefficients}{plsone}(object, ...)
\method{residuals}{plsone}(object, ...)
\method{fitted}{plsone}(object, ...)
\method{formula}{plsone}(x, ...)
}
\arguments{
  \item{x}{an object of class inheriting from 'plsone'}
  \item{object}{an object of class inheriting from 'plsone'}
  \item{formula}{an object of class \code{formula} (or one that can be coerced to that class):
  a symbolic description of the model to be fitted. The details of model specification are given under Details.}
  \item{data}{an optional data frame, list or environment (or object coercible by as.data.frame to a data frame)
  containing the variables in the model. If not found in data, the variables are taken from environment(formula)}
  \item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}
  \item{na.action}{a function which indicates what should happen when the data contain NAs. The default is set
  by the na.action setting of options, and is na.fail if that is unset. The 'factory-fresh' default is na.omit.
  Another possible value is NULL, no action. Value na.exclude can be useful. Value na.pass is used
  to take in account the missing values.}
  \item{y}{a numeric vector corresponding to the response}
  \item{X}{a numeric data.frame corresponding to explanatory variables (see \code{model.matrix})}
  \item{nf}{an integer indicating the number of kept components}
  \item{scale}{a logical value. If scale is TRUE, X is scaled by dividing each variable
  by its sample standard deviation (by deafult, scale=TRUE).}
  \item{method}{the method to be used in fitting the model. The default method "plsone.fit" uses partial least squares (PLS).
  The only current alternative is "model.frame" which returns the model frame and does no fitting.}
  \item{type}{the type of coefficient. The default "components" retruns the coefficients of regression of y on the components.
  The alternative "varaibles" gives coefficients based on initial explanatory variables}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  The function 'plsone' is regression methods based on partail least square.\cr
  The default method 'plsone.fit' uses NIPALS algorithm.\cr
  The Coordinates ch,th, ph and ch correspond to Slopes equation in linear least squares regression without intercept
  term and can be compute with missing data.\cr
  + additional details on the computation of the terms.\cr
  + The library propose two functuion to fit model: plsone.fit programmed only in R and cplsone.fit based on C programs.\cr
  The function plsone is base on the function delete.intercept from the library 'pls'.\cr
}
\value{
  The function 'plsone' returns an object of class inheriting from 'plsone'.
  \item{lm}{an object 'lm' based on the PLS components (th)}
  \item{bh}{an numeric vector which contains coefficients associated with the explanatory variables}
  \item{y}{a numeric vector corresponding to the response}
  \item{x}{a numeric data.frame corresponding to explanatory variables}
  \item{ch}{a numeric vector containing regression coefficients of y on the components 'th'}
  \item{th}{an object 'data.frame' containing the conponents 'th'.}
  \item{uh}{an object 'data.frame' containing normalised residuals for h > 1.}
  \item{ph}{an object 'data.frame' containing regression coefficient  of the \eqn{h_{th}} columns of \eqn{X_{h-1}} on t(h).}
  \item{wh}{an object 'data.frame' containing regression coefficient  of the \eqn{h_{th}} columns of \eqn{X_{h-1}} on y(h-1).}
  \item{whx}{an object 'data.frame' containing regression coefficent of th on explanatory variables.}
  \item{Q2}{a numeric vector conatining the criteria Q2 given by: \eqn{Q2=1-PRESS_{h}/RSS_{h-1}}.}
  \item{nf}{an integer indicating the number of kept components.}
  \item{formula}{the initial formula supplied.}
  \item{method}{the name of the fitter function used, currently always "plsone.fit". The C version is called "cplsone.fit".}
  \item{call}{generally match.call()}
}
\references{
Wold H. (1966) Estimation of principal components and related models by iterative least squares.
In P. Krishnaiah, editors.Multivariate Analysis, Academic Press, 391–420.\cr
Wold H. (1975) Modeling in Complex Situations with Soft Information, Third World
Congress of Econometric Society, August 21-26, Toronto, Canada.\cr
Wold S., Martens H. & Wold H. (1983) « The multivariate calibration problem in chemistry
solved by the PLS method», in  Proc. Conf. Matrix Pencils, Ruhe, A. & Kågstrøm, B. (Eds),
March 1982, Lecture Notes in Mathematics, Springer Verlag, Heidelberg, pp. 286- 293.\cr
Tenenhaus M.(1998) La Regression PLS. Theorie et pratique. Technip, Paris.\cr
}
\note{
The function 'pls' use the function 'delete.intercept' from the package \code{\link[pls:pls-package]{pls}}.
}
\seealso{
\code{\link[pls:plsr]{plsr}}, \code{\link[pls:delete.intercept]{delete.intercept}},
\code{\link[pls:mvr]{mvr}}, \code{\link[stats:lm]{lm}}}
\examples{
require(pls)
data(yarn)
yarn.pls <- plsr(density ~ NIR, 6, data = yarn, validation = "CV")
plstest1 <- plsone(density ~ NIR, nf=6, data = yarn,scale=FALSE)
cor(plstest1$th,yarn.pls$scores)
cor(fitted(plstest1),yarn.pls$fitted.values[,,6])
plot(fitted(plstest1),yarn.pls$fitted.values[,,6],pch=20,panel.first=c(grid()),cex=2)
abline(0,1,col="red",lwd=2)
# update fonction
plstest2 <- update(plstest1,nf=2)
plstest2
plstest1
# with missing data
data(yarn)
yarn$NIR[sample(1:(268*28),300)] <- rep(NA,300)
plsna1 <- plsone(density ~ NIR, nf=6, data = yarn,scale=FALSE,na.action=na.pass)
plot(fitted(plstest1),fitted(plsna1),pch=20,panel.first=c(grid()),cex=2)
abline(0,1,col="red",lwd=2)
}


\name{gplsone}
\alias{gplsone}
\alias{gplsone.fit}
\alias{cgplsone.fit}
\alias{residuals.gplsone}
\alias{scalew}
\title{GPLS model based on NIPALS algorithm}
\description{
Generalized Partial least squares regresssion based on NIPALS algorithm (IWPLS).
}
\usage{
gplsone(formula, data,subset,weights=NULL,na.action,offset,
  family=gaussian,nf = 2,method="gplsone.fit",scale=FALSE,
  tol=1e-9,deps=1e-20,nitermax=100, ...)
gplsone.fit(y,X,nf=2,family,weights=NULL,tol=1e-9,deps=1e-20,
  nitermax=100,scale=FALSE,...)
cgplsone.fit(y,X,nf=2,family,weights=NULL,tol=1e-9,deps=1e-20,
  nitermax=100,scale=FALSE,...)
\method{residuals}{gplsone}(object,...)
}
\arguments{
  \item{object}{an object of class inheriting from 'gplsone'.}
  \item{formula}{an object of class \code{formula} (or one that can be coerced to that class):
  a symbolic description of the model to be fitted. The details of model specification are given under Details.}
  \item{data}{an optional data frame, list or environment (or object coercible by as.data.frame to a data frame)
  containing the variables in the model. If not found in data, the variables are taken from environment(formula)}
  \item{family}{a description of the error distribution and link function to be used in the model.
  This can be a character   string naming a family function, a family function or the result
  of a call to a family function. (See "family" for details of family functions.)}
  \item{weights}{an optional vector of weights to be used in the fitting process. Should be NULL or a numeric vector.}
  \item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}
  \item{na.action}{a function which indicates what should happen when the data contain NAs. The default is set
  by the na.action setting of options, and is na.fail if that is unset. The 'factory-fresh' default is na.omit.
  Another possible value is NULL, no action. Value na.exclude can be useful. Value na.pass is used
  to take in account the missing values.}
  \item{offset}{this can be used to specify an a priori known component to be included in the linear predictor
  during fitting. This should be NULL or a numeric vector of length either one or equal to the number of cases.
  One or more offset terms can be included in the formula instead or as well, and if both are specified
  their sum is used (See model.offset).}
  \item{y}{a numeric vector corresponding to the response}
  \item{X}{a numeric data.frame corresponding to explanatory variables (see \code{model.matrix})}
  \item{nf}{an integer indicating the number of kept components}
  \item{scale}{a logical value. If scale is TRUE, X is scaled by dividing each variable
  by its sample standard deviation (by deafult, scale=TRUE).}
  \item{method}{the method to be used in fitting the model. The default method "plsone.fit" uses partial least squares (PLS).
  The only current alternative is "model.frame" which returns the model frame and does no fitting.}
 \item{tol}{positive convergence tolerance (by default, tol=1e-9).}
 \item{deps}{positive tolerance (used in the computation of convergence value,by default deps=1e-20).}
 \item{nitermax}{integer giving the maximal number of IWPLS iterations (by default, nitermax=100).}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  The function 'plsone' is regression methods based on partail least square.\cr
  The default method 'plsone.fit' uses NIPALS algorithm.\cr
  The Coordinates ch,th, ph and ch correspond to Slopes equation in linear least squares regression without intercept
  term and can be compute with missing data.\cr
  + additional details on the computation of the terms.\cr
  The function 'gplsone' is base on two functions: scalew (from the library ade4) and delete.intercept (from the library pls).
  For more details on these both methods, you can directly consult their documentations.\cr
  The functions 'cgplsone.fit' is based on C code and the function 'gpls' is programmed in R.
}
\value{
  The function 'gplsone' returns an object of class inheriting from 'plsone'.
  \item{glm}{an object 'glm' based on the PLS components (th)}
  \item{bh}{an numeric vector which contains coefficients associated with the explanatory variables}
  \item{y}{a numeric vector corresponding to the response}
  \item{x}{a numeric data.frame corresponding to explanatory variables}


  \item{ch}{a numeric vector containing regression coefficients of y on the components 'th'}
  \item{th}{an object 'data.frame' containing the conponents 'th'.}
  \item{uh}{an object 'data.frame' containing normalised residuals for h > 1.}
  \item{ph}{an object 'data.frame' containing regression coefficient  of the \eqn{h_{th}} columns of \eqn{X_{h-1}} on t(h).}
  \item{wh}{an object 'data.frame' containing regression coefficient  of the \eqn{h_{th}} columns of \eqn{X_{h-1}} on y(h-1).}
  \item{whx}{an object 'data.frame' containing regression coefficent of th on explanatory variables.}

???? \cr

  \item{nf}{an integer indicating the number of kept components.}
  \item{formula}{the initial formula supplied.}
  \item{method}{the name of the fitter function used, currently always "plsone.fit".}
  \item{call}{generally match.call()}
}
\references{
Marx 1996
}
\note{
The function 'pls' use the function 'delete.intercept' from the package \code{\link[pls:pls-package]{pls}}.
}
\seealso{
\code{\link[gpls:gpls]{gpls}}, \code{\link[pls:delete.intercept]{delete.intercept}},
\code{\link[stats:glm]{glm}}}
\examples{
require(pls)
data(yarn)

}


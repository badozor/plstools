\name{predict.plsone}
\alias{predict.plsone}
\title{Predict method for Partial Least Squares Fits}
\description{
The function 'predict.plsone'retruns predicted values,
confidence and tolerance intervals.
}
\usage{
\method{predict}{plsone}(object, newdata, se.fit = FALSE,interval="none",level = 0.95,...)
}
\arguments{
  \item{object}{Object of class inheriting from 'plsone'}
  \item{newdata}{An optional data frame in which to look for variables with which to predict.
  If omitted, the fitted values are used.}
  \item{se.fit}{A switch indicating if standard errors are required.}
  \item{interval}{none, confidence, prediction}
  \item{level}{Tolerance/confidence level}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
add formula for tolerance and confidence interval
}
\value{
The function 'predict.plsone'retruns ...
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{
Tenenhaus M. (1998) La Regression PLS. Theorie et pratique. Technip, Paris.\cr}
\seealso{\code{\link{plsone}}, \code{\link[stats:predict.lm]{predict.lm}}}
\examples{
require(pls)
data(yarn)
yarn.pls <- plsr(density ~ NIR, 6, data = yarn, validation = "CV")
plstest1 <- plsone(density ~ NIR, nf=6, data = yarn,scale=FALSE)
plot(predict(plstest1)$fitted,plstest1$y)
}


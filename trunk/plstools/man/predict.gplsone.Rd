\name{predict.gplsone}
\alias{predict.gplsone}
\title{Predict method for Generalized Partial Least Squares Fits}
\description{
The function 'predict.gplsone'retruns predicted values,
confidence and tolerance intervals.
}
\usage{
\method{predict}{gplsone}(object,newdata=NULL,type="link",se.fit=FALSE,interval="none",level=0.95,...)
}
\arguments{
  \item{object}{Object of class inheriting from 'plsone'}
  \item{newdata}{An optional data frame in which to look for variables with which to predict.
  If omitted, the fitted values are used.}
  \item{type}{the type of prediction required. The default is on the scale of the linear predictors;
  the alternative "response" is on the scale of the response variable. Thus for a default binomial model
  the default predictions are of log-odds (probabilities on logit scale) and type = "response" gives
  the predicted probabilities. The "terms" option returns a matrix giving the fitted values of each term
  in the model formula on the linear predictor scale. The value of this argument can be abbreviated.}
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
\seealso{\code{\link{gplsone}}, \code{\link[stats:predict.lm]{predict.glm}}}
\examples{
require(pls)
data(yarn)

}


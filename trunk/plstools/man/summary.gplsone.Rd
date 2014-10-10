\name{summary.gplsone}
\alias{summary.gplsone}
\alias{print.summary.gplsone}
\title{Summarizing Generalized Partial Least Squares model}
\description{
This function provides additional elements to evaluate the goodness of fit
and the quality of the model (object of class inheriting from 'gplsone').
}
\usage{
\method{summary}{gplsone}(object, ...)
\method{print}{summary.gplsone}(x, ...)
}
\arguments{
  \item{object}{an object of class inheriting from 'gplsone'}
  \item{x}{an object of class inheriting from 'summary.gplsone'}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
 decrire le calcul des elements ...\cr

+ additional details on the computation of each elemnts\cr
+ At the moment this function is very slow to compute the diagnostic elements.\cr
\cr
}
\value{
The function 'summary.gplsone' retruns an list of class 'summary.plsone' containing
several elements associated with the model diagnostic:\cr

  \item{bh}{an numeric vector which contains coefficients associated with the explanatory variables}
  \item{n}{observation number}
  \item{p}{number of explanatory varaibles.}
  \item{nf}{an integer indicating the number of kept components (h).}
}
\references{
Marx 1996 \cr}
\seealso{\code{\link{gplsone}},\code{\link[stats:summary.glm]{summary.glm}}}
\examples{
require(pls)
data(yarn)

}



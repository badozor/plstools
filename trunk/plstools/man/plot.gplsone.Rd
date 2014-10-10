\name{plot.gplsone}
\alias{plot.gplsone}
\title{Plot Diagnostics for an gplsone Object}
\description{

in construction

}
\usage{
\method{plot}{gplsone}(x,xax=1,yax=2,...)
}
\arguments{
  \item{x}{an object of class inheriting from 'gplsone'}
  \item{xax}{the column number for the x-axis}
  \item{yax}{the column number for the y-axis}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
The selected plots are drawn on a graphics device.\cr
1: Observed vs Expected values\cr
2: Normal QQ-plot of residuals\cr
3: Representation of Components th (in 2 dimensions, by default xax=1 and yax=2).\cr
4: Prediction vs the \eqn{h^{th}} component (by default xax=1).\cr
6: Representation of Weights (in 2 dimensions, by default xax=1 and yax=2).\cr
9: Representation of correlations with the conponents
(in 2 dimensions, by default xax=1 and yax=2) \cr
}
\value{
  x is invisibly returned.
}
\references{
Tenenhaus M.(1998) La Regression PLS. Theorie et pratique. Technip, Paris.\cr}
\seealso{ \code{\link{gplsone}}, \code{\link[stats:qqnorm]{qqnorm}},
\code{\link[ade4:s.arrow]{s.arrow}},
\code{\link[ade4:s.label]{s.label}}, \code{\link[stats:plot.glm]{plot.glm}}}
\examples{
require(pls)
data(yarn)

}


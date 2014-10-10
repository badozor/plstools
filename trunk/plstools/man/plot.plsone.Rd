\name{plot.plsone}
\alias{plot.plsone}
\title{Plot Diagnostics for an plsone Object}
\description{
Nine plots (selectable by which.plot) are currently available:
Observed vs Expected values, Normal QQ-plot of residuals,
Representation of components th, Prediction vs the \eqn{h^{th}} component,
Representation of DModXN, Representation of Weights,
Representation of DModYN, Representationof Hotellin T2,
Representation of correlations with the conponents.
}
\usage{
\method{plot}{plsone}(x, xax = 1, yax = 2, mfrow = NULL, which.plot =1:9,...)
}
\arguments{
  \item{x}{an object of class inheriting from 'plsone'}
  \item{xax}{the column number for the x-axis}
  \item{yax}{the column number for the y-axis}
  \item{mfrow}{a vector of the form 'c(nr,nc)', otherwise computed by as special
  own function n2mfrow}
  \item{which.plot}{a numeric vector containing the numbers of The selected plots
   (see details)}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
The selected plots are drawn on a graphics device.\cr
1: Observed vs Expected values\cr
2: Normal QQ-plot of residuals\cr
3: Representation of Components th (in 2 dimensions, by default xax=1 and yax=2).\cr
4: Prediction vs the \eqn{h^{th}} component (by default xax=1).\cr
5: Representation of DModXN\cr
6: Representation of Weights (in 2 dimensions, by default xax=1 and yax=2).\cr
7: representation of DModYN\cr
8: Representationof Hotellin T2\cr
9: Representation of correlations with the conponents
(in 2 dimensions, by default xax=1 and yax=2) \cr
}
\value{
  x is invisibly returned.
}
\references{
Tenenhaus M.(1998) La Regression PLS. Theorie et pratique. Technip, Paris.\cr}
\seealso{ \code{\link{plsone}}, \code{\link[stats:qqnorm]{qqnorm}},
\code{\link[ade4:s.arrow]{s.arrow}},
\code{\link[ade4:s.label]{s.label}}, \code{\link[stats:plot.lm]{plot.lm}}}
\examples{
require(pls)
data(yarn)
plstest <- plsone(density ~ NIR, nf=6, data = yarn,scale=FALSE)
plot(plstest)
plot(plstest,which.plot=1:2)
plot(plstest,which.plot=c(1:2,6))
}


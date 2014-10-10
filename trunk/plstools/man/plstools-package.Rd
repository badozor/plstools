\name{plstools-package}
\alias{plstools-package}
\alias{plstools}
\docType{package}
\title{
PLS regression based on NIPALS algoritm
}
\description{
This package contains PLS regression functions based on NIPALS algoritm
and some tools to evaluate the goodness of fit and the model quality.
}
\details{
\tabular{ll}{
Package: \tab plstools\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2008-09-18\cr
License: \tab GPL version 2 or newer\cr
LazyLoad: \tab yes\cr
Depends: \tab R (>= 2.7.1),ade4\cr
Suggests: \tab ade4, boot\cr
}
The function 'pls' use the function 'delete.intercept' from the package \code{\link[pls:pls-package]{pls}}.
}
\author{
Pierre Bady <pierre.bady@unil.ch>\cr
Maintainer: Pierre Bady <pierre.bady@unil.ch>\cr
}
\references{
Wold H. (1966) Estimation of principal components and related models by iterative least squares.
In P. Krishnaiah, editors.Multivariate Analysis, Academic Press, 391–420.\cr
Wold H. (1975) Modeling in Complex Situations with Soft Information, Third World
Congress of Econometric Society, August 21-26, Toronto, Canada.\cr
Wold S., Martens H. & Wold H. (1983) « The multivariate calibration problem in chemistry
solved by the PLS method», in  Proc. Conf. Matrix Pencils, Ruhe, A. & Kågstrøm, B. (Eds),
March 1982, Lecture Notes in Mathematics, Springer Verlag, Heidelberg, pp. 286- 293.\cr
Tenenhaus M.(1998) La Regression PLS. Theorie et pratique. Technip, Paris.\cr}
\seealso{
\code{\link[pls:pls-package]{pls}},
\code{\link[gpls:gpls-package]{gpls}},
\code{\link[ade4:ade4-package]{ade4}},
\code{\link[boot:boot-package]{boot}}
}
\examples{
require(pls)
data(yarn)
yarn.pls <- plsr(density ~ NIR, 6, data = yarn, validation = "CV")
plstest1 <- plsone(density ~ NIR, nf=6, data = yarn,scale=FALSE)
cor(plstest1$th,yarn.pls$scores)
cor(fitted(plstest1),yarn.pls$fitted.values[,,6])
plot(fitted(plstest1),yarn.pls$fitted.values[,,6],pch=20,panel.first=c(grid()),cex=2)
abline(0,1,col="red")
abline(0,1,col="red",lwd=2)
plstest2 <- update(plstest1,nf=2)
plstest2
plstest1
summary(plstest1)
summary(plstest2)
plot(plstest1,which.par=1:2)
plot(plstest1,which.par=c(1:2,6))
}

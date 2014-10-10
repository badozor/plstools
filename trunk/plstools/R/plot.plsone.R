`plot.plsone` <-
function (x, xax = 1, yax = 2, mfrow = NULL, which.plot =1:9,...)
{
    if (!inherits(x, "plsone"))
      stop("non convenient argument!")
    lm1 <- x$lm
    s1 <- summary(x)
    opar <- par(ask = par("ask"), mfrow = par("mfrow"))
    on.exit(par(opar))
    if (is.null(mfrow))
        mfrow <- n2mfrow(length(which.plot))
    par(mfrow = mfrow)
    if (length(which.plot) > prod(mfrow))
        par(ask = TRUE)
  plotFittedY <- function(fitted,y,...){
    plot(fitted, y, xlab = "expected values", ylab = "observed values",
        pch = 17,panel.first=c(grid()))
    abline(lm(y ~ fitted), lwd = 2, col = "red")
    abline(0, 1, lty = 2, lwd = 1.5)
    title("observed and expected values")
    }
  plotQQresid <- function(lm1,...){
    qqnorm(residuals(lm1, type = "response"), pch = 17,panel.first=c(grid()))
    qqline(residuals(lm1, type = "response"), lwd = 2, col = "red")
    abline(0, 1, lty = 2, lwd = 1.5)
    }
  plotComponents <- function(th,xax,yax,...)
    s.label(th, xax = xax, yax = yax, sub = "components Th")
  plotHotellingT2 <- function(s1,...){
    n <- s1$n
    seuil <- s1$FT2[s1$nf]
    Tw <- s1$T2[, s1$nf]
    max1 <- max(c(seuil, Tw))
    max1 <- max1 + 0.1 * max1
    plot(Tw, ylab = "T2", type = "h", ylim = c(0, max1),panel.first=c(grid()))
    abline(h = seuil, lty = 3, col = "red")
    title("Hotelling T2")
    w <- Tw[Tw >= seuil]
    if (length(w) != 0) {
        badvalues <- as.data.frame(cbind(c(1:n)[Tw >= seuil],
            Tw[Tw >= seuil] + 0.1 * max(Tw)))
        text(badvalues, label = badvalues[, 1])
        }
    }
  plotDModXN <- function(s1,...){
    n <- s1$n
    seuil <- s1$DcritX
    nw <- s1$DModXN
    max1 <- max(c(seuil, nw))
    max1 <- max1 + 0.1 * max1
    plot(nw, ylab = "DModXN", type = "h", ylim = c(0, max1),panel.first=c(grid()))
    abline(h = seuil, lty = 3, col = "red")
    title("DModXn")
    w <- nw[nw >= seuil]
    if (length(w) != 0) {
        badvalues <- as.data.frame(cbind(c(1:n)[nw >= seuil],
            nw[nw >= seuil] + 0.1 * max(nw)))
        text(badvalues, label = badvalues[, 1])
        }
    }
  plotDModYN<- function(s1,...){
    n <- s1$n
    seuil <- s1$DcritY
    nw <- s1$DModYN
    max1 <- max(c(seuil, nw))
    max1 <- max1 + 0.1 * max1
    plot(nw, ylab = "DModYN", type = "h", ylim = c(0, max1),panel.first=c(grid()))
    abline(h = seuil, lty = 3, col = "red")
    title("DModYn")
    w <- nw[nw >= seuil]
    if (length(w) != 0) {
        badvalues <- as.data.frame(cbind(c(1:n)[nw >= seuil],
            nw[nw >= seuil] + 0.1 * max(nw)))
        text(badvalues, label = badvalues[, 1])
      }
    }
  plotCorrComp <- function(rxt2,xax,yax,...)
    s.arrow(rxt2, xax = xax, yax = yax, sub = "Correlation with components")
  plotWeights <- function(whx,xax,yax,...)
    s.arrow(whx, xax = xax, yax = yax, sub = "Weights")
  plotFittedComp <- function(x,xax,...){
    fit <- predict(x,se.fit=TRUE,interval="confidence")
    pred1 <- fit$fitted
    ec1 <- fit$lower
    ec2 <- fit$upper
    coox <- x$th[, xax]
    ord1 <- order(coox)
    plot(coox[ord1], pred1[ord1], type = "n", ylim = c(min(ec1),
        max(ec2)), ylab = "expected values", xlab = paste("Th[,", xax, "]",
        sep = ""),panel.first=c(grid()))
    polygon(cbind(c(coox[ord1], rev(coox[ord1])), c(ec1[ord1],
        rev(ec2[ord1]))), border = FALSE, col = gray(0.5))
    lines(coox[ord1], pred1[ord1], lwd = 1)
    title("component vs prediction")
    }
    for (i in which.plot) {
        switch(i,
        plotFittedY(lm1$fitted,x$y),
        plotQQresid(lm1),
        plotComponents(x$th,xax,yax),
        plotFittedComp(x,xax),
        plotDModXN(s1),
        plotWeights(x$whx,xax,yax),
        plotDModYN(s1),
        plotHotellingT2(s1),
        plotCorrComp(s1$rxt2,xax,yax))
    }
    invisible(match.call())
}

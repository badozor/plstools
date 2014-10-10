`plsone.fit` <-
function (y, X, nf = 2, scale=TRUE, ...)
{
    if (!is.data.frame(X))
        stop("data.frame expected")
    fip <- function(X, y, k) {
        x <- X[k, ]
        X <- X[-k, ]
        y <- y[-k]
        nr1 <- nrow(X)
        nc1 <- ncol(X)
        to <- rep(0, nr1)
        w <- rep(1/sqrt(nc1), nc1)
        for (j in 1:nc1) {
            ye <- y[!is.na(X[, j])]
            w[j] <- sum(X[, j] * y, na.rm = TRUE)/sum(ye * ye, na.rm = TRUE)
        }
        w <- w/sqrt(sum(w * w, na.rm = TRUE))
        for (i in 1:nr1) {
            we <- w[!is.na(X[i, ])]
            to[i] <- sum(X[i, ] * w, na.rm = TRUE)/sum(we * we, na.rm = TRUE)
        }
        co <- sum(y * to)/sum(to * to)
        return(co * sum(x * w, na.rm = TRUE))
    }
    Xh <- scale(X, center=TRUE,scale=scale)
    yh <- ys <- as.vector(scale(y))
    nc <- ncol(X)
    nr <- nrow(X)
    Ch <- Wh <- Uh <- NULL
    Th <- Ph <- Q2 <- NULL
    yest <- rep(0, nr)
    for (h in 1:nf) {
        rssh <- sum(yh * yh)
        pressh <- 0
        for (r in 1:nr) {
            pressh <- pressh + (yh[r] - fip(Xh, yh, r))^2
        }
        Q2 <- c(Q2, 1 - pressh/rssh)
        th <- rep(0, nr)
        ph <- rep(0, nc)
        wh <- rep(1/sqrt(nc), nc)
        for (j in 1:nc) {
            yhj <- yh[!is.na(Xh[, j])]
            wh[j] <- sum(Xh[, j] * yh, na.rm = TRUE)/sum(yhj * yhj, na.rm = TRUE)
        }
        wh <- wh/sqrt(sum(wh * wh, na.rm = TRUE))
        for (i in 1:nr) {
            whj <- wh[!is.na(Xh[i, ])]
            th[i] <- sum(Xh[i, ] * wh, na.rm = TRUE)/sum(whj * whj, na.rm = TRUE)
        }
        for (j in 1:nc) {
            ph[j] <- sum(Xh[, j] * th, na.rm = TRUE)/sum(th * th, na.rm = TRUE)
        }
        Xh <- Xh - th %*% t(ph)
        ch <- sum(yh * th, na.rm = TRUE)/sum(th * th, na.rm = TRUE)
        uh <- yh/ch
        yh <- yh - ch * th
        Ch <- c(Ch, ch)
        Wh <- cbind(Wh, wh)
        Th <- cbind(Th, th)
        Ph <- cbind(Ph, ph)
        Uh <- cbind(Uh, uh)
    }
    Whx <- Wh %*% solve(crossprod(Ph, Wh))
    bh <- Whx %*% Ch
    res <- list()
    Th <- as.data.frame(Th)
    names(Th) <- paste("t", 1:nf, sep = "")
    row.names(Th) <- row.names(X)
    Ph <- as.data.frame(Ph)
    names(Ph) <- paste("cs", 1:nf, sep = "")
    row.names(Ph) <- names(X)
    Wh <- as.data.frame(Wh)
    names(Wh) <- paste("cs", 1:nf, sep = "")
    row.names(Wh) <- names(X)
    lm1 <- lm(y ~ ., data = Th)
    mx <- apply(X, 2, function(x) mean(x, na.rm = TRUE))
    if(scale)
      sdx <- apply(X, 2, function(x) sd(x, na.rm = TRUE))
    else
      sdx <- rep(1,ncol(X))
    my <- mean(y, na.rm = TRUE)
    sdy <- sd(y, na.rm = TRUE)
    a <- sum(bh * (-mx)/sdx) * sdy + my
    b <- bh * sdy/sdx
    coefx <- c(a, b)
    names(coefx) <- c("(Intercept)", names(X))
    res <- list(lm = lm1, bh = coefx, y = y, x = X, ch = Ch,
        th = Th, uh = Uh, ph = Ph, wh = Wh, whx = Whx, Q2 = Q2, 
        nf = nf, call = match.call())
    class(res) <- "plsone"
    return(res)
}

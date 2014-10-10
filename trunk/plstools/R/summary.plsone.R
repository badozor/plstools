`summary.plsone` <- function (object, ...) 
{
    if (!inherits(object, "plsone"))
        stop("object 'plsone' expected")
    n <- length(object$y)
    nf <- object$nf
    th <- object$th
    df <- as.data.frame(scale(object$x))
    y <- scale(object$y)
    p <- ncol(df)
    w <- as.data.frame(na.omit(cbind(df, th)))
    df1 <- w[, c(1:p)]
    th1 <- w[, -c(1:p)]
    res <- list()
    rx2 <- ry2 <- NULL
    r2 <- function(x1, y1) (cor(x1, y1))^2
    for (i in 1:nf) {
        tw <- th1[, i]
        w <- unlist(lapply(df1, function(x) r2(x, tw)))
        rx2 <- cbind(rx2, w)
        ry2 <- c(ry2, r2(y, th[, i]))
    }
    res$rx2 <- rx2
    res$ry2 <- ry2
    sh2 <- unlist(lapply(object$th, var))
    T2 <- NULL
    for (i in 1:nf) {
        w <- sapply(1:n, function(x) sum((th[x, 1:i]^2)/sh2[1:i]) * 
            n/(n - 1))
        T2 <- cbind(T2, w)
    }
    res$T2 <- T2
    Ts2 <- n * (n - nf) * T2/(nf * (n * n - 1))
    res$Ts2 <- Ts2
    f1 <- function(nf, n) nf * (n * n - 1) * qf(0.95, nf, n - 
        nf)/(n * (n - nf))
    res$FT2 <- unclass(sapply(1:nf, function(x) f1(x, n)))
    rxt2 <- ryt2 <- NULL
    for (i in 1:nf) {
        tw <- th1[, i]
        w <- unlist(lapply(df1, function(x) cor(x, tw)))
        rxt2 <- cbind(rxt2, w)
        ryt2 <- c(ryt2, cor(y, th[, i]))
    }
    res$rxt2 <- rxt2
    res$ryt2 <- ryt2
    fei <- function(x, th) return((x - predict(lm(x ~ ., data = th), 
        newdata = th))^2)
    w <- as.data.frame(lapply(df, function(x) fei(x, th)))
    si <- sqrt(apply(w, 1, function(x) sum(x, na.rm = TRUE))/(p -
        nf)) * sqrt(n/(n - nf - 1))
    s0 <- sqrt(sum(w, na.rm = TRUE)/((p - nf) * (n - nf - 1)))
    res$DModX <- matrix(si)
    res$DModXN <- matrix(si/s0)
    fi <- lm(y ~ ., data = th)$residuals^2
    fi <- abs(sqrt(fi) * sqrt(n/(n - nf - 1)))
    res$DModY <- fi
    res$DModYN <- fi/sqrt(sum(fi^2)/(n - nf - 1))
    fout <- function(x, coef = 1.5) {
        stats <- stats::fivenum(x, na.rm = TRUE)
        iqr <- diff(stats[c(2, 4)])
        out <- (stats[4] + coef * iqr)
        return(out)
    }
    res$DcritX <- fout(res$DModXN, coef = 1.5)
    res$DcritY <- fout(res$DModYN, coef = 1.5)
    s1 <- summary(object$lm)
    res$coefficients <- s1$coefficients
    res$r.squared <- s1$r.squared
    res$adj.r.squared <- s1$adj.r.squared
    res$bh <- object$bh
    res$n <- n
    res$nf <- nf
    res$p <- p
    res$call <- object$call
    class(res) <- "summary.plsone"
    return(res)
}

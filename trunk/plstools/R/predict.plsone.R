`predict.plsone` <-
function(object, newdata, se.fit = FALSE,interval="none",level = 0.95,...)
{
    if (!inherits(object, "plsone"))
        stop("object 'plsone' expected")
    n <- nrow(object$th)
    nf <- object$nf
    ddl <- n - nf - 1
    rss <- sqrt(sum(object$lm$residuals^2)/object$lm$df.residual)
#########
    if (missing(newdata)) {
        w <- apply(data.frame(lapply(object$th, function(x) (x * 
            x)/sum(x * x))), 1, sum)
        if(interval=="prediction")
          sd2 <- sqrt(w + 1/n + 1)
        else
          sd2 <- sqrt(w + 1/n)
        sefit <- rss * sd2
        fitted <- object$lm$fitted
    }
    else {
###### à revoir !!!
#    mt <- attr(mf, "terms")
#    df <- as.data.frame(delete.intercept(model.matrix(mt, mf)))


#########
        df <- newdata
        bh <- object$bh
        fitted <- as.vector(as.matrix(df) %*% as.matrix(bh[-1]) + 
            bh[1])
        mx <- unlist(lapply(object$x, function(x) mean(x, na.rm = TRUE)))
        sdx <- unlist(lapply(object$x, function(x) sd(x, na.rm = TRUE)))
        tw <- NULL
        for (i in 1:nrow(df)) {
            w <- unlist(apply(object$whx, 2, function(x) sum(x * 
                (df[i, ] - mx)/sdx, na.rm = TRUE)))
            tw <- rbind(tw, w)
        }
        w <- sapply(1:nf, function(x) (tw[, x]^2)/sum(object$th[, 
            x]^2))
        w <- apply(w, 1, sum)
## revoir pour cette option ?
        if(interval=="prediction")
          sd2 <- sqrt(w + 1/n + 1)
        else
          sd2 <- sqrt(w + 1/n)
        sefit <- rss * sd2
        names(fitted) <- names(df)
    }
    res <- list()
    res$fitted <- fitted
    if(se.fit)
      res$se.fit <- sefit
    if(interval=="confidence" | interval=="prediction"){
      res$lower <- fitted - sefit*qt(1 - (1 - level)/2, df = ddl)
      res$upper <- fitted + sefit*qt(1 - (1 - level)/2, df = ddl)
      }
    return(res)
}

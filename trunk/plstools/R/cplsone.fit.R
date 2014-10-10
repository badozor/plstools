`cplsone.fit` <-
function (y, X, nf = 2, scale=TRUE,...)
{
  Xh <- scale(X, center=TRUE,scale=scale)
  yh <- as.vector(scale(y))
  p <- ncol(X)
  n <- nrow(X)
  Ch <- rep(0,nf)
  Wh <- matrix(0,nrow=p,ncol=nf)
  Ph <- matrix(0,nrow=p,ncol=nf)
  Th <- matrix(0,nrow=n,ncol=nf)
  Uh <- matrix(0,nrow=n,ncol=nf)
  fit1 <- .C("plsFitter",as.integer(p),as.integer(n),as.integer(nf),as.double(Xh),as.double(yh),
    as.double(Ch),as.double(Wh),as.double(Th),as.double(Ph),as.double(Uh),NAOK = TRUE,PACKAGE="plstools")
  Ch <- fit1[[6]]
  Wh <- matrix(fit1[[7]],nrow=p,ncol=nf)
  Ph <- matrix(fit1[[9]],nrow=p,ncol=nf)
  Th <- matrix(fit1[[8]],nrow=n,ncol=nf)
  Uh <- matrix(fit1[[10]],nrow=n,ncol=nf)
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
  if(scale){
    sdx <- apply(X, 2, function(x) sd(x, na.rm = TRUE))
  }else{ sdx <- rep(1,ncol(X))}
  my <- mean(y, na.rm = TRUE)
  sdy <- sd(y, na.rm = TRUE)
  a <- sum(bh * (-mx)/sdx) * sdy + my
  b <- bh * sdy/sdx
  coefx <- c(a, b)
  names(coefx) <- c("(Intercept)", names(X))
  res <- list(lm = lm1, bh = coefx, y = y, x = X, ch = Ch,
        th = Th, uh = Uh, ph = Ph, wh = Wh, whx = Whx,
        nf = nf, call = match.call())
  class(res) <- "plsone"
  return(res)
}


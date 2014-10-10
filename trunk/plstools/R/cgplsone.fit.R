`cgplsone.fit` <-
function(y,X,nf=2,family,weights=NULL,tol=1e-9,deps=1e-20,nitermax=100,scale=FALSE,...){
  if (!is.data.frame(X))
    stop("data.frame expected")
# fonctions utilitaires
  varmu<- family$variance
  dlink <- family$mu.eta
  link <- family$linkfun
  linkinv <- family$linkinv
# pour la gestion des NA =>  à faire
  nc <- ncol(X)
  nobs <- nr <- nrow(X)
  if(is.null(weights))
    weights <-  rep(1, nr)
# initialisation à change en fonction de la loi !!! 10/18/08 12:58:19
  mustart <- NULL
etastart <- NULL
  eval(family$initialize)
  fk <- link(mustart)
#  print("fk:")
#  print(fk)
  mu <- linkinv(fk)
  v <- (dlink(fk)^2/varmu(mu))
#  print("v:")
#  print(v)
  Ek <- scalew(X,weights=as.numeric(v),center=TRUE,scale=scale)
# conteneurs
  Tk <- matrix(0,nrow=nr,ncol=nf)
  Wk <- matrix(0,nrow=nc,ncol=nf)
  Qk <- rep(0,nf)
  Pk <- matrix(0,nrow=nc,ncol=nf)
  Uk <- matrix(0,nrow=nr,ncol=nf)
#
  beta <- rep(1,nc)
  betaold <- beta/1000
  dbeta <- 1
  niter <- 1
  convergence <- FALSE
  while( (dbeta > tol)&& (niter < nitermax)){
    fit1 <- .C("plsFitterwt",as.integer(nc),as.integer(nr),as.integer(nf),as.double(Ek),as.double(fk),
      as.double(v),as.double(Qk),as.double(Wk),as.double(Tk),as.double(Pk),as.double(Uk),NAOK = TRUE,PACKAGE="plstools")
    Qk <- fit1[[7]]
    Wk <- matrix(fit1[[8]],nrow=nc,ncol=nf)
    Pk <- matrix(fit1[[10]],nrow=nc,ncol=nf)
    Tk <- matrix(fit1[[9]],nrow=nr,ncol=nf)
    Uk <- matrix(fit1[[11]],nrow=nr,ncol=nf)
    eta <- as.numeric(weighted.mean(fk,v)+Tk%*%Qk)
    #
    mu <- linkinv(eta)
    v <- (dlink(eta)^2/varmu(mu))
    if(sum(v=="NaN")>0)
      break
    fk <- eta + (1/dlink(eta))*(y-mu)
    Ek <- scalew(X,weights=as.numeric(v) ,center=TRUE,scale=scale)
    niter <- niter+1
    betaold <- beta
    Whx <- Wk %*% solve(crossprod(Pk, Wk))
    beta <- Whx %*% Qk
    dbeta <- max(abs(beta - betaold)/(abs(betaold) + deps))
    }
  if(dbeta < tol)
    convergence <- TRUE
  glm1 <- glm(y~as.matrix(Tk),family=family)
# coeff en fonction des variables
    Whx <- Wk %*% solve(crossprod(Pk, Wk))
    bh <- Whx %*% Qk
    mx <- apply(X, 2, function(x) mean(x, na.rm = TRUE))
  meta <- mean(eta,na.rm = TRUE)
  a <- meta - sum(mx*bh)
  coefx <- c(a,bh)
  names(coefx) <- c("(Intercept)",names(X))
#results
  res <- list(glm = glm1,bh = coefx,y = y,x = X,qh = Qk,th = as.data.frame(Tk),
    ph = Pk,wh = Wk,nf = nf,niter = niter,convergence = convergence,call = match.call())
  #res$uk <- Uk
  class(res) <- c("gplsone","plsone")
  return(res)
}

`gplsone.fit` <-
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
  mu <- linkinv(fk)
  v <- (dlink(fk)^2/varmu(mu))
  Ek <- scalew(X,weights=as.numeric(v),center=TRUE,scale=scale)
# conteneurs
  Tk <- matrix(0,nrow=nr,ncol=nf)
  Wk <- matrix(0,nrow=nc,ncol=nf)
  Qk <- rep(0,nf)
  Pk <- matrix(0,nrow=nc,ncol=nf)
  #Uk <- matrix(0,nrow=nr,ncol=nf)
  beta <- rep(1,nc)
  betaold <- beta/1000
  dbeta <- 1
  niter <- 1
  convergence <- FALSE
  while( (dbeta > tol)&& (niter < nitermax)){
# debut iteration pour les tk
    for(k in 1:nf){
        tk <- rep(0,nr)
        pk <- rep(0,nc)
        wk <- rep(1/sqrt(nc),nc)
# début algorithme nipals
        for(j in 1:nc){
          fkj <- fk[!is.na(Ek[,j])]
          vj <- v[!is.na(Ek[,j])]
          wk[j] <- sum(Ek[,j]*fk*v,na.rm=TRUE)/sum(fkj*vj*fkj,na.rm=TRUE)
          }
        for(i in 1:nr){
          wkj <- wk[!is.na(Ek[i,])]
          tk[i] <- sum(Ek[i,]*wk,na.rm=TRUE)/sum(wkj*wkj,na.rm=TRUE)
          }
         for(j in 1:nc){
          pk[j] <- sum(Ek[,j]*v*tk,na.rm=TRUE)/sum(tk*v*tk,na.rm=TRUE)
         }
        wk <- wk/sqrt(sum(wk*wk,na.rm=TRUE))
        qk <- sum(fk*tk*v,na.rm=TRUE)/sum(tk*v*tk,na.rm=TRUE)
        fk <- fk - qk*tk
        Ek <- Ek - tk%*%t(pk)
        Tk[,k] <- tk
        Qk[k] <- qk
        Pk[,k] <- pk
        Wk[,k] <- wk
        #Uk[,k] <- fk/qk
# fin algorithme nipals
        }
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

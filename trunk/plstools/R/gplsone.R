`gplsone` <-
function (formula, data,subset,weights=NULL,na.action,offset,family=gaussian,
  nf = 2,method="gplsone.fit",scale=FALSE,tol=1e-9,deps=1e-20,nitermax=100, ...)
{
    call <- match.call()
    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
      }
    if (missing(data))
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action","offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    switch(method, model.frame = return(mf), gplsone.fit = 1,cgplsone.fit = 2, stop("invalid 'method': ", method))
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm))
            names(Y) <- nm
      }
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(Y), 0)
    X <- as.data.frame(delete.intercept(X))
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights))
        stop("'weights' must be a numeric vector")
    offset <- as.vector(model.offset(mf))
    if (!is.null(weights) && any(weights < 0))
        stop("negative weights not allowed")
    if (!is.null(offset)) {
        if (length(offset) == 1)
            offset <- rep(offset, NROW(Y))
        else if (length(offset) != NROW(Y))
            stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                length(offset), NROW(Y)), domain = NA)
        }
    fit <- gplsone.fit(Y, X, nf = nf, family=family,weights=weights,tol=tol,deps=deps,
                nitermax=nitermax,scale=scale,...)
                
    if(method=="gplsone.fit"){
      fit <- gplsone.fit(Y, X, nf = nf, family=family,weights=weights,tol=tol,deps=deps,
                nitermax=nitermax,scale=scale,...)
    }else{
      fit <- cgplsone.fit(Y, X, nf = nf, family=family,weights=weights,tol=tol,deps=deps,
                nitermax=nitermax,scale=scale,...)
    }
    fit$formula <- formula
    fit$method <- method
    fit$call <- call
    return(fit)
}

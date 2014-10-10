`plsone` <-
function (formula, data,subset,na.action,nf = 2,scale=TRUE,method="plsone.fit",...)
{
    call <- match.call()
    if (missing(data))
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset","na.action"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    switch(method, model.frame = return(mf), plsone.fit = 1,cplsone.fit = 2, stop("invalid 'method': ", method))
    mt <- attr(mf, "terms")
    ry <- model.response(mf, "numeric")
    rdf <- as.data.frame(delete.intercept(model.matrix(mt, mf)))
    if(method=="plsone.fit"){
      fit <- plsone.fit(ry, rdf, nf = nf, scale, ...)
    }else{
      fit <- cplsone.fit(ry, rdf, nf = nf, scale, ...)
    }
    fit$formula <- formula
    fit$method <- method
    fit$call <- call
    return(fit)
}

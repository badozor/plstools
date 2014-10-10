`scalew` <-
function (X, weights = rep(1, nrow(X)), center = TRUE, scale = TRUE, na.rm=TRUE)
{
# d'apres scaleweigths de D.CHESSEL issu du package ade4
    X <- as.matrix(X)
    n <- nrow(X)
    if (length(weights) != n)
        stop("length of weights must equal the number of rows in x")
    if (any(weights < 0) || (s <- sum(weights)) == 0)
        stop("weights must be non-negative and not all zero")
    weights <- weights/s
    center <- if (center)
        apply(weights * X, 2, function(x) sum(x,na.rm=na.rm))
    else 0
    X <- sweep(X, 2, center)
    norm <- apply(X * X * weights, 2, function(x) sum(x,na.rm=na.rm))
    norm[norm <= 1e-07 * max(norm)] <- 1
    if (scale)
        X <- sweep(X, 2, sqrt(norm), "/")
    return(X)
}

`print.plsone` <-
function(x,...){
    if (!inherits(x, "plsone"))
        stop("object 'plsone' expected")
    if (!is.null(x$lm))
      nomAnalysis <- "\nPLS model (first generation)\n"
    if (!is.null(x$glm))
      nomAnalysis <- "\nGPLS model (first generation)\n"
    cat(nomAnalysis)
    cat("\n$call: ")
    cat(class(x))
    cat("\n$class: ")
    print(x$call)
    cat("$nf:", x$nf, "components used\n")
    cat("\nformula with components:\n\n")
    if (!is.null(x$lm))
        coefw <- x$lm$coefficients
    if (!is.null(x$glm))
        coefw <- x$glm$coefficients
    print(coefw)
    cat("\n")
}

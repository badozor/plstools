`print.summary.plsone` <-
function (x,digits = max(3, getOption("digits") - 3), ...)
{
    if (!inherits(x, "summary.plsone"))
        stop("x 'plsone' expected")
    cat("Partial Least Square Regression\n")
    print(x$call)
    cat("coefficients based on components:\n")
    print(round(x$coefficients,digits))
    cat("Multiple R-squared: ",round(x$r.squared,digits),"\t")
    cat("Adjusted R-squared: ",round(x$adj.r.squared,digits),"\n")
    cat("threshold values :\n")
    cat("DcritX: ", round(x$DcritX,digits), "\n")
    cat("DcritY: ", round(x$DcritY,digits), "\n")
    cat("F(Hotelling's T2):\n")
    print(round(x$FT2,digits))
    cat("\n")
    invisible(x)
}


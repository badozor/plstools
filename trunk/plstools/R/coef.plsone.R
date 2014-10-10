`coef.plsone` <-
function(object,type="components",...){
  if (!inherits(object,"plsone"))
    stop("object 'plsone' expected")
  if(type=="components"){
    if(!is.null(object$glm))
      w <- object$glm$coefficients
    else w <- object$lm$coefficients
    }
  else if(type=="variables")
    w <- object$bh
  else
    stop("non convenient 'type'!")
  return(w)
}

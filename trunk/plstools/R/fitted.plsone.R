`fitted.plsone` <-
function(object,...){
  if(!is.null(object$glm))
    w <- object$glm$fitted
  else w <- object$lm$fitted
  w
}

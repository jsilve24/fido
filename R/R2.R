##' Generic Method to Calculate R2 for Fitted Model
##'
##' @param m model object
##' @param ... other arguments to pass
##' @return vector
##' @name r2
NULL

##' @rdname r2
##' @export
r2 <- function(m, ...) {
 UseMethod("r2", m) 
}

r2_internal <- function(eta.hat, eta){
  resid <- eta - eta.hat
  ss.res <- apply(resid, c(3), function(x) sum(diag(var(x))))
  eta.mean <- apply(eta, c(1,3), mean)
  ss.tot <- apply(eta, c(3), function(x) sum(diag(var(x))))

  return(1-ss.res/ss.tot)
}


##' @rdname r2
##' @details
##'   Calculates Posterior over Linear Model R2 as:
##'   \deqn{1-\frac{SS_{res}}{SS_tot}}
##'   where \eqn{SS} is defined in terms of trace of variances
##' @export
r2.pibblefit <- function(m, ...){
  req(m, c("Lambda", "Eta", "X"))

  eta.hat <- predict(m, newdata=NULL, response="LambdaX")
  r2_internal(eta.hat, m$Eta)
}


##' @rdname r2
##' @details
##'   Calculates Posterior over Basset Model R2 as:
##'   \deqn{1-\frac{SS_{res}}{SS_tot}}
##' @export
r2.bassetfit <- function(m, ...){
  req(m, c("Lambda", "Eta"))

  eta.hat <- predict(m, newdata=NULL, response="Lambda")
  r2_internal(eta.hat, m$Eta)
}

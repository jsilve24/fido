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
  ## variance / sum of squares being defined as total variation
  ## (trace of covariance matrix)
  var.res <- apply(resid, c(3), function(x) sum(diag(var(t(x)))))
  var.tot <- apply(eta, c(3), function(x) sum(diag(var(t(x)))))

  ## Note as we are doing this all in the latent space (under the multinomial)
  ## it guarantees that 1-var.res/var.tot is equal to var_fit/(var_fit+var_res)
  ## exactly as discussed in Gelman 2019
  return(1-var.res/var.tot)
}


##' @rdname r2
##' @details
##'   Calculates Posterior over Linear Model R2 as:
##'   \deqn{1-\frac{SS_{res}}{SS_{tot}}}
##'   where \eqn{SS} is defined in terms of trace of variances
##' 
##'   Method of calculating R2 is multivariate version of the Bayesian R2 proposed
##'   by Gelman, Goodrich, Gabry, and Vehtari, 2019
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
##' 
##'   Method of calculating R2 is multivariate version of the Bayesian R2 proposed
##'   by Gelman, Goodrich, Gabry, and Vehtari, 2019
##' @export
r2.bassetfit <- function(m, ...){
  req(m, c("Lambda", "Eta"))

  eta.hat <- predict(m, newdata=NULL, response="Lambda")
  r2_internal(eta.hat, m$Eta)
}

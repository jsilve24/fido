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
  var.pred <- apply(eta.hat, c(3), function(x) sum(diag(var(t(x)))))

  return(var.pred/(var.pred + var.res))
}

##' @rdname r2
##' @param covariates vector of indices for covariates to include in calculation of R2 (default:NULL
##'   means include all covariates by default). When non-null, all covariates not specified are set
##'   to zero for prediction.
##' @details Calculates Posterior over Linear Model R2 as: \deqn{1-\frac{SS_{res}}{SS_{tot}}} where
##'   \eqn{SS} is defined in terms of trace of variances
##' 
##'   Method of calculating R2 is multivariate version of the Bayesian R2 proposed
##'   by Gelman, Goodrich, Gabry, and Vehtari, 2019
##' @export
r2.pibblefit <- function(m, covariates=NULL, ...){
  req(m, c("Lambda", "Eta", "X"))
  newdata <- m$X
  if (!is.null(covariates)){
    ## Defensive
    covariates <- unique(covariates)
    stopifnot("covariates must be integer valued" = all(round(covariates)==covariates))
    stopifnot("some passed covariates outside of range 1:Q" = max(covariates)<=m$Q & min(covariates)>=1)
    ## END DEFENSE

    ## set non-included covariates to zero

    covariates.exclude <- setdiff(1:m$Q, covariates)
    newdata[covariates.exclude,] <- 0
  }
  eta.hat <- predict(m, newdata=newdata, response="LambdaX")
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

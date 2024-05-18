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

row_var <- function(x) {
        rowSums((x - rowMeans(x))^2) / (dim(x)[2] - 1)
}

r2_internal <- function(eta.hat, eta) {
        resid <- eta - eta.hat
        ## variance / sum of squares being defined as total variation
        ## (trace of covariance matrix)
        var.res <- apply(resid, c(3), function(x) sum(row_var(x)))
        var.pred <- apply(eta.hat, c(3), function(x) sum(row_var(x)))

        return(var.pred / (var.pred + var.res))
}


##' @rdname r2
##' @param covariates vector of indices for covariates to include in calculation of R2 (default:NULL
##'   means include all covariates by default). When non-null, all covariates not specified are set
##'   to zero for prediction.
##' @details Calculates Posterior over Linear Model R2 as: \deqn{1-\frac{SS_{res}}{SS_{tot}}}{1-(SS_(res)/SS_(tot))} where
##'   \eqn{SS} is defined in terms of trace of variances
##'
##'   Method of calculating R2 is multivariate version of the Bayesian R2 proposed
##'   by Gelman, Goodrich, Gabry, and Vehtari, 2019
##' @export
r2.pibblefit <- function(m, covariates = NULL, ...) {
        req(m, c("Lambda", "Eta", "X"))
        newdata <- m$X
        if (!is.null(covariates)) {
                ## Defensive
                covariates <- unique(covariates)
                stopifnot("covariates must be integer valued" = all(round(covariates) == covariates))
                stopifnot("some passed covariates outside of range 1:Q" = max(covariates) <= m$Q & min(covariates) >= 1)
                ## END DEFENSE

                ## set non-included covariates to zero

                covariates.exclude <- setdiff(1:m$Q, covariates)
                newdata[covariates.exclude, ] <- 0
        }
        eta.hat <- predict(m, newdata = newdata, response = "LambdaX")
        r2_internal(eta.hat, m$Eta)
}


##' @rdname r2
##' @param covariates vector of indices for covariates to include in calculation of R2 (default:NULL
##'   means include all covariates by default). When non-null, all covariates not specified are set
##'   to zero for prediction.
##' @param components vector of indices for components of the GP model to include in the calculation of R2, i.e. which
##' elements in the list of Theta/Gamma should be used for calculating R2 (default:NULL
##' means to include all components by default). When non-null, all components not specified are removed
##' for prediction.
##' @details
##'   Calculates Posterior over Basset Model R2 as:
##'   \deqn{1-\frac{SS_{res}}{SS_{tot}}}{1-(SS_(res)/SS_(tot))}
##'
##'   Method of calculating R2 is multivariate version of the Bayesian R2 proposed
##'   by Gelman, Goodrich, Gabry, and Vehtari, 2019
##' @export
r2.bassetfit <- function(m, covariates = NULL, components = NULL, ...) {
    req(m, c("Lambda", "Eta"))
    newdata <- m$X
    m.used <- m
    if (!is.null(covariates)) {
      ## Defensive
      
      ## Finding the linear component
      linear.comp <- which(sapply(m$Gamma, is.matrix))
      stopifnot("there must be a linear component to use the covariates input" =  !is.null(linear.comp))
      Q <- dim(m$Lambda[[linear.comp]])[2]
      
      covariates <- unique(covariates)
      stopifnot("covariates must be integer valued" = all(round(covariates) == covariates))
      stopifnot("some passed covariates outside of range 1:Q" = max(covariates) <= Q & min(covariates) >= 1)
      ## END DEFENSE
      
      covariates.exclude <- setdiff(1:Q, covariates)
      linear.comp <- which(sapply(m$Gamma, is.matrix))
      m.used$Lambda[[linear.comp]][,covariates.exclude,] <- 0
    }
    
    if(!is.null(components)){
      components <- unique(components)
      stopifnot("components must be integer valued" = all(round(components) == components))
      stopifnot("some passed components outside of range 1:length(m$Lambda)" = max(components) <= length(m$Lambda) & min(components) >= 1)
      ## END DEFENSE
      # components.exclude <- setdiff(1:length(m$Lambda), components)
      # for(i in components.exclude){
      #   m.used$Lambda[[i]] <- array(0, dim(m$Lambda[[i]]))
      # }
      m.used$Lambda <- m.used$Lambda[components]
      m.used$Gamma <- m$Gamma[components]
      m.used$Theta <- m$Theta[components]
    }
    eta.hat <- predict(m.used, newdata = newdata, response = "Lambda")
    r2_internal(eta.hat, m$Eta)
}

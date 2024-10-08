# internal function 
new_maltipoofit <- function(D, N, Q, P, coord_system, iter=NULL, 
                            alr_base=NULL, ilr_base=NULL,
                            Eta=NULL, Lambda=NULL,Sigma=NULL, Sigma_default=NULL, 
                            Y=NULL, X=NULL, upsilon=NULL, 
                            Theta=NULL, Xi=NULL,Xi_default=NULL, Gamma=NULL, 
                            init=NULL, ellinit=NULL, names_categories=NULL, names_samples=NULL, 
                            names_covariates=NULL, VCScale=NULL, U=NULL){
  m <- new_pibblefit(D, N, Q, coord_system, iter, alr_base, ilr_base,
                     Eta, Lambda, Sigma, Sigma_default, 
                     Y, X, upsilon, Theta, Xi,Xi_default, Gamma, 
                     init, ellinit, names_categories, names_samples)
  m$VCScale <- VCScale
  m$U <- U
  m$ellinit <- ellinit
  m$P <- P
  m$names_covariates <- names_covariates
  class(m) <- c("maltipoofit", "pibblefit")
}


#' Create maltipoofit object 
#' 
#' @inheritParams pibblefit
#' @inheritParams maltipoo_fit
#' @param VCScale scale factors (delta) for variance components 
#' @param P number of variance components
#' @return object of class maltipoofit
#' @seealso \code{\link{maltipoo}}
#' @noRd
maltipoofit <- function(D, N, Q, P, coord_system, iter=NULL,  
                        alr_base=NULL, ilr_base=NULL,
                        Eta=NULL, Lambda=NULL, Sigma=NULL, Sigma_default=NULL, 
                        Y=NULL, X=NULL, upsilon=NULL, 
                        Theta=NULL, Xi=NULL,Xi_default=NULL, Gamma=NULL, 
                        init=NULL, ellinit=NULL, names_categories=NULL, names_samples=NULL, 
                        names_covariates=NULL, VCScale=NULL, U=NULL){
  m <- new_maltipoofit(D, N, Q, coord_system, iter, alr_base, ilr_base,
                    Eta, Lambda, Sigma, Sigma_default, 
                    Y, X, upsilon, Theta, Xi,Xi_default, Gamma, 
                    init, ellinit, names_categories, names_samples, 
                    names_covariates, VCScale, U)
  verify_maltipoofit(m)
  return(m)
}



#' Simple verification of passed multipoo object
#' Changed from \code{verify.maltipoofit} to \code{verify_maltipoofit}.
#' This is due to an issue with Generics and S3 methods. 
#' Before, R CMD check result in the NOTE (CRAN doesn't like most NOTES):
#' "Apparent methods for exported generics not registered"
#' @param m an object of class multipoo
#' @param ... not used
#' @return throws error if any verification tests fail
#' @noRd
verify_maltipoofit <- function(m,...){
  verify.pibblefit(m)
  stopifnot(is.integer(m$P))
  ifnotnull(m$VCScale, check_dims(m$VCScale, m$P, "VCScale"))
  ifnotnull(m$U, check_dims(m$U, c(m$P*m$Q, m$Q), "U"))
  ifnotnull(m$ellinit, check_dims(m$ellinit, m$P, "ellinit"))
}

#' require elements to be non-null in pibblefit or throw error
#' Changed from \code{req.maltipoofit} to \code{req_maltipoofit}.
#' This is due to an issue with Generics and S3 methods. 
#' Before, R CMD check result in the NOTE (CRAN doesn't like most NOTES):
#' "Apparent methods for exported generics not registered"
#' 
#' @inheritParams req
#' @return Throws an error if null
#' @noRd
req_maltipoofit <- function(m, r){
  present <- sapply(m[r], is.null)
  if(any(present)){
    stop("maltipoofit object does not contain required components:", r[present])
  }
}
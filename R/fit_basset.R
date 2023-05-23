
#' Interface to fit basset models
#' 
#' Basset (A Lazy Learner) - non-linear regression models in fido
#'
#' @param Y D x N matrix of counts (if NULL uses priors only)
#' @param X Q x N matrix of covariates (cannot be NULL)
#' @param upsilon dof for inverse wishart prior (numeric must be > D) 
#'   (default: D+3)
#' @param Theta A function from dimensions dim(X) -> (D-1)xN (prior mean of gaussian process). For an additive GP model, can be a list of functions from dimensions dim(X) -> (D-1)xN + a (optional) matrix of size (D-1)xQ for the prior of a linear component if desired.
#' @param Gamma A function from dimension dim(X) -> NxN (kernel matrix of gaussian process). For an additive GP model, can be a list of functions from dimension dim(X) -> NxN + a QxQ prior covariance matrix if a linear component is specified. It is assumed that the order matches the order of Theta.
#' @param Xi (D-1)x(D-1) prior covariance matrix
#'   (default: ALR transform of diag(1)*(upsilon-D)/2 - this is 
#'   essentially iid on "base scale" using Aitchison terminology)
#' @param init (D-1) x Q initialization for Eta for optimization 
#' @param pars character vector of posterior parameters to return
#' @param m object of class bassetfit 
#' @param ... other arguments passed to \link{pibble} (which is used internally to 
#'  fit the basset model)
#' 
#' @details the full model is given by:
#'    \deqn{Y_j \sim Multinomial(\Pi_j)} 
#'    \deqn{\Pi_j = \Phi^{-1}(\Eta_j)}
#'    \deqn{\Eta \sim MN_{D-1 x N}(\Lambda, \Sigma, I_N)}
#'    \deqn{\Lambda \sim GP_{D-1 x Q}(\Theta(X), \Sigma, \Gamma(X))}
#'    \deqn{\Sigma \sim InvWish(\upsilon, \Xi)}
#'  Where Gamma(X) is short hand for the Gram matrix of the Kernel function. 
#'  
#'  Alternatively can be used to fit an additive GP of the form:
#'    \deqn{Y_j \sim Multinomial(\Pi_j)} 
#'    \deqn{\Pi_j = \Phi^{-1}(\Eta_j)}
#'    \deqn{\Eta \sim MN_{D-1 x N}(\Lambda, \Sigma, I_N)}
#'    \deqn{\Lambda = \Lambda_1 + ... + \Lambda_p + \Beta X}
#'    \deqn{\Lambda_1 \sim GP_{D-1 x Q}(\Theta_1(X), \Sigma, \Gamma_p(X))}
#'    ...
#'    \deqn{\Lambda_p \sim GP_{D-1 x Q}(\Theta_1(X), \Sigma, \Gamma_1(X))}
#'    \deqn{\Beta \sim MN(\Theta_B, \Sigma, \Gamma_B)}
#'    \deqn{\Sigma \sim InvWish(\upsilon, \Xi)}
#'  Where Gamma(X) is short hand for the Gram matrix of the Kernel function. 
#'  
#'  Default behavior is to use MAP estimate for uncollaping the LTP 
#'  model if laplace approximation is not preformed. 
#' @return an object of class bassetfit
#' @md
#' @name basset_fit
NULL

#' @rdname basset_fit
#' @export
basset <- function(Y=NULL, X, upsilon=NULL, Theta=NULL, Gamma=NULL, Xi=NULL, 
                   init=NULL, pars=c("Eta", "Lambda", "Sigma"), ...){
  
  args <- list(...)
  ncores <- args_null("ncores", args, -1)
  seed <- args_null("seed", args, sample(1:2^15, 1))
  ret_mean <- args_null("ret_mean", args, FALSE)
  n_samples <- args_null("n_samples", args, 2000)
  
  D <- nrow(Y)
  if (is.null(upsilon)) upsilon <- D+3  # default is minimal information 
  # but with defined mean
  if (is.null(Xi)) {
    # default is iid on base scale
    # G <- cbind(diag(D-1), -1) ## alr log-constrast matrix
    # Xi <- 0.5*G%*%diag(D)%*%t(G) ## default is iid on base scale
    Xi <- matrix(0.5, D-1, D-1) # same as commented out above 2 lines
    diag(Xi) <- 1               # same as commented out above 2 lines
    Xi <- Xi*(upsilon-D) # make inverse wishart mean Xi as in previous lines 
  }
  
  
  ## adding functionality so that Theta and Gamma can be a list
  if(typeof(Theta) == "list" | typeof(Gamma) == "list"){
    if(typeof(Gamma) != "list" | typeof(Theta) != "list"){
      stop("Theta and Gamma must both be lists if one of them is a list.")
    }
    if(length(Gamma) != length(Theta)){
      stop("Theta and Gamma must be of the same length.")
    }
    
    ## theta
    Theta_trans <- list()
    for(i in 1:length(Theta)){
      if(is.matrix(Theta[[i]])){
        Theta_trans[[i]] <- Theta[[i]] %*% X
      } else if(is.function(Theta[[i]])){
        Theta_trans[[i]] <- Theta[[i]](X)
      } else{
      stop(paste("Element ", i, "of Theta not supported!"))
    }}
    Theta_comb <- Reduce('+', Theta_trans)
    ## gamma
    Gamma_trans <- list()
    for(i in 1:length(Gamma)){
      if(is.matrix(Gamma[[i]])){
        Gamma_trans[[i]] <- t(X) %*% Gamma[[i]] %*% X
      } else if(is.function(Gamma[[i]])){
        Gamma_trans[[i]] <- Gamma[[i]](X)
      } else{
        stop(paste("Element ", i, "of Gamma not supported!"))
      }
    }
    
    Gamma_comb <- Reduce('+', Gamma_trans)
    
    ## fitting the joint model
    collapse_samps <- pibble(Y, X=diag(ncol(X)), upsilon, Theta_comb, Gamma_comb, Xi, init, pars, ...)
    
    ## fitting uncollapse using the joint samples
    
    Lambda <- list()
    
    ## Sampling for each of the additive components
    l <- length(Theta_trans)
    for(i in 1:l){
      if(i == 1){
        samp_mean <- Reduce("+", Theta_trans[-c((i+1):l)])
        eta_samples <- sweep(collapse_samps$Lambda, c(1,2), samp_mean, FUN = "-")
      } else{
        red_lambda <- Reduce("+", Lambda)
        samp_mean <- array(NA, dim = dim(red_lambda))
        for(j in 1:dim(samp_mean)[3]){
          samp_mean[,,j] <- red_lambda[,,j] + Reduce("+", Theta_trans[-c((i+1):l)])
        }
        eta_samples <- collapse_samps$Lambda - samp_mean
      }
      Gamma_comb_red <- Reduce('+', Gamma_trans[c(i:l)])
      
      if(i != l){
        if(is.matrix(Theta[[i]])){
          fitu <- uncollapsePibble_sigmaKnown(eta_samples, X, Theta[[i]], Gamma[[i]], Gamma_comb_red, Xi, collapse_samps$Sigma, upsilon, 
                                   ret_mean=ret_mean, linear = TRUE, ncores=ncores, seed=seed)
          LambdaX <- array(NA, dim = c(nrow(fitu$Lambda), ncol(X), n_samples))
          for(j in 1:n_samples){
            LambdaX[,,j] <- fitu$Lambda[,,j] %*% X
          }
          Lambda[[i]] <- LambdaX
        } else{
          fitu <- uncollapsePibble_sigmaKnown(eta_samples, diag(ncol(X)), Theta_trans[[i]], Gamma_trans[[i]], Gamma_comb_red, Xi, collapse_samps$Sigma, upsilon, 
                                                     ret_mean=ret_mean, linear = FALSE, ncores=ncores, seed=seed)
          Lambda[[i]] <- fitu$Lambda
        }
      } else{
        row.names(collapse_samps$Lambda) <- NULL
        Lambda[[i]] <- collapse_samps$Lambda - samp_mean
      }
    }
    ## pretty output ##
    out <- list()
    if ("Eta" %in% pars){
      out[["Eta"]] <- collapse_samps$Eta
    }
    if ("Lambda" %in% pars){
      out[["Lambda"]] <- Lambda
    }
    if ("Sigma" %in% pars){
      out[["Sigma"]] <- fitu$Sigma
    }
    
    # By default just returns all other parameters
    out$N <- collapse_samps$N
    out$Q <- nrow(X)
    out$D <- collapse_samps$D
    out$Y <- Y
    out$X <- X
    out$upsilon <- upsilon
    out$Theta <- Theta
    out$Xi <- Xi
    out$Gamma <- Gamma
    out$init <- collapse_samps$init
    out$iter <- dim(collapse_samps$Eta)[3]
    # for other methods
    out$names_categories <- rownames(Y)
    out$names_samples <- colnames(Y)
    out$names_covariates <- rownames(X)
    out$coord_system <- "alr"
    out$alr_base <- collapse_samps$D
    out$summary <- NULL
    out$logMarginalLikelihood <- collapse_samps$logMarginalLikelihood
    attr(out, "class") <- c("pibblefit")
    # add names if present 
    # if (use_names) out <- name(out)

  } else{
    if (!is.null(Theta)) {
      Theta_train <- Theta(X)
    } else {
      Theta <- function(X) matrix(0, nrow(Y)-1, ncol(X)) 
      Theta_train <- Theta(X)
    }
    if (!is.null(Gamma)) {
      Gamma_train <- Gamma(X)
    } else {
      stop("No Default Kernel For Gamma Implemented")
    }
    out <- pibble(Y, X=diag(ncol(X)), upsilon, Theta_train, Gamma_train, Xi, init, pars, ...)
    out$X <- X
    out$Q <- as.integer(nrow(X))
  }
  out$Theta <- Theta
  out$Gamma <- Gamma
  class(out) <- c("bassetfit", "pibblefit")
  verify(out)
  return(out)
}

#' @rdname basset_fit
#' @export
refit.bassetfit <- function(m, pars=c("Eta", "Lambda", "Sigma"), ...){
  # Store coordinates and tranfsorm to cannonical representation
  l <- store_coord(m)
  m <- to_alr(m, m$D)
  
  # Concatenate parameters to pass to basset function
  argl <- list(...)
  argl$pars <- pars
  ml <- as.list(m)
  argl <- c(ml, argl)
  
  # Need to handle iter as part of m but no n_samples passed
  # in this situation should pull iter from m and pass as n_samples to pibble 
  if (is.null(argl[["n_samples"]]) & !is.null(m$iter)) argl[["n_samples"]] <- m$iter 
  if (typeof(m$Theta) == "list"){argl$Q <- ncol(m$X)}
  # pass to basset function
  m <- do.call(basset, argl)
  
  # Reapply original coordinates
  m <- reapply_coord(m, l)
  verify(m)
  return(m)
}

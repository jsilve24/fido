
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
#' @param linear A vector denoting which rows of X should be used if  a linear component is specified. Default is all rows.
#' @param init (D-1) x Q initialization for Eta for optimization 
#' @param pars character vector of posterior parameters to return
#' @param m object of class bassetfit 
#' @param newdata Default is `NULL`. If non-null, newdata is used in the uncollapse sampler in place of X.
#' @param ... other arguments passed to \link{pibble} (which is used internally to 
#'  fit the basset model)
#' 
#' @details the full model is given by:
#'    \deqn{Y_j \sim Multinomial(\pi_j)}{Y_j \sim Multinomial(Pi_j)}  
#'    \deqn{\pi_j = \Phi^{-1}(\eta_j)}{Pi_j = Phi^(-1)(Eta_j)}
#'    \deqn{\eta \sim MN_{D-1 \times N}(\Lambda, \Sigma, I_N)}{Eta \sim MN_(D-1 x N)(Lambda, Sigma, I_N)}
#'    \deqn{\Lambda \sim GP_{D-1 \times Q}(\Theta(X), \Sigma, \Gamma(X))}{Lambda \sim GP_{D-1 times Q}(Theta(X), Sigma, Gamma(X))}
#'    \deqn{\Sigma \sim InvWish(\upsilon, \Xi)}{Sigma \sim InvWish(upsilon, Xi)}
#'  Where \eqn{\Gamma(X)}{Gamma(X)}  is short hand for the Gram matrix of the Kernel function. 
#'  
#'  Alternatively can be used to fit an additive GP of the form:
#'    \deqn{Y_j \sim Multinomial(\pi_j)}{Y_j \sim Multinomial(Pi_j)} 
#'    \deqn{\pi_j = \Phi^{-1}(\eta_j)}{Pi_j = Phi^(-1)(Eta_j)}
#'    \deqn{\eta \sim MN_{D-1 \times N}(\Lambda, \Sigma, I_N)}{Eta \sim MN_(D-1 times N)(Lambda, Sigma, I_N)}
#'    \deqn{\Lambda = \Lambda_1 + ... + \Lambda_p + B X}{Lambda = Lambda_1 + ... + Lambda_p + Beta X}
#'    \deqn{\Lambda_1 \sim GP_{D-1 \times Q}(\Theta_1(X), \Sigma, \Gamma_1(X))}{Lambda_1 \sim GP_{D-1 x Q}(Theta_1(X), Sigma, Gamma_1(X))}
#'    \deqn{...}
#'    \deqn{\Lambda_p \sim GP_{D-1 \times Q}(\Theta_p(X), \Sigma, \Gamma_p(X))}{Lambda_p \sim GP_{D-1 x Q}(Theta_p(X), Sigma, Gamma_p(X))}
#'    \deqn{B \sim MN(\Theta_B, \Sigma, \Gamma_B)}{Beta \sim MN(Theta_B, Sigma, Gamma_B)}
#'    \deqn{\Sigma \sim InvWish(\upsilon, \Xi)}{Sigma \sim InvWish(upsilon, Xi)}
#'  Where \eqn{\Gamma(X)}{Gamma(X)} is short hand for the Gram matrix of the Kernel function. 
#'  
#'  Default behavior is to use MAP estimate for uncollaping the LTP 
#'  model if laplace approximation is not preformed. 
#' @return an object of class bassetfit
#' @md
#' @name basset_fit
NULL

#' @rdname basset_fit
#' @export
basset <- function(Y=NULL, X, upsilon=NULL, Theta=NULL, Gamma=NULL, Xi=NULL, linear = NULL,
                   init=NULL, pars=c("Eta", "Lambda", "Sigma"), newdata = NULL, ...){
  
  args <- list(...)
  ncores <- args_null("ncores", args, -1)
  seed <- args_null("seed", args, sample(1:2^15, 1))
  ret_mean <- args_null("ret_mean", args, FALSE)
  n_samples <- args_null("n_samples", args, 2000)

  D <- nrow(Y)
  N <- ncol(Y)
  if(ncol(X) != N) stop("The number of columns in X and Y must match.")
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
  
  if(is.null(linear)){
    message("No rows of X were specified. Using all rows...")
    linear <- 1:nrow(X)
  } else{
    ## Check that linear doesn't contain values less than 1 or greater than nrow(X)
    if(min(linear) < 1 | max(linear) > nrow(X)){
      stop("Please verify that all values of 'linear' are between 1 and nrow(X).")
    }
  }
  
  
  ## adding functionality so that Theta and Gamma can be a list
  if(typeof(Theta) == "list" | typeof(Gamma) == "list"){
    if(typeof(Gamma) != "list" | typeof(Theta) != "list"){
      stop("Theta and Gamma must both be lists if one of them is a list.")
    }
    if(length(Gamma) != length(Theta)){
      stop("Theta and Gamma must be of the same length.")
    }
    
    ## evaluating theta and gamma
    ## theta
    theta_eval <- function(Theta, X, linear){
      if(is.matrix(Theta)){
        Q.red <- length(linear)
        if(ncol(Theta) != Q.red) stop("The dimension of the linear Theta element is incorrect! Please ensure it matches the dimensions of the desired linear components.")
        return(Theta %*% X[linear,])
      } else if(is.function(Theta)){
        return(Theta(X))
      } else{
        stop("An element of Theta is not supported! All elements must be a matrix or list.")
      }
    }
    
    ## gamma
    gamma_eval <- function(Gamma, X, linear){
      if(is.matrix(Gamma)){
        Q.red <- length(linear)
        if(length(linear) == 1){
          X.red <- matrix(X[linear,], nrow =1)
        } else{
          X.red <- X[linear,]
        }
        if(ncol(Gamma) != Q.red | nrow(Gamma) != Q.red) stop("The dimension of the linear component of Gamma element is incorrect! Please ensure it matches the dimensions of the desired linear components.")
        return(t(X.red) %*% Gamma %*% X.red)
      } else if(is.function(Gamma)){
        return(Gamma(X))
      } else{
        stop("An element of Gamma is not supported! All elements must be a matrix or list.")
      }
    }
    Theta_trans <- list()
    Gamma_trans <- list()
    for(i in 1:length(Theta)){
      Theta_trans[[i]] <- theta_eval(Theta[[i]], X, linear)
      Gamma_trans[[i]] <- gamma_eval(Gamma[[i]], X, linear)
    }
    
    Theta_comb <- Reduce('+', Theta_trans)
    Gamma_comb <- Reduce('+', Gamma_trans)
    
    ## fitting the joint model
    ## newdata auto handled by pibble
    collapse_samps <- pibble(Y, X=diag(ncol(X)), upsilon, Theta_comb, Gamma_comb, Xi, init, pars, newdata = newdata, ...)
    
    ## fitting uncollapse using the joint samples
    Lambda <- list()
    Lambda.out <- list()
    
    ##Updating the number of samples... Useful if it returns the MAP estimate
    if(dim(collapse_samps$Eta)[3] != n_samples){
      warning("Using MAP estimates for uncollapsing...")
    } 
    
    ## Sampling for each of the additive components
    num.comp <- length(Theta_trans)
    Lambda <- list()
    Lambda.out <- list()
    
    add.uncollapse <- function(coll_samples, X, Theta, Gamma, Gamma_comb, Xi, Sigma, upsilon, ret_mean, ncores, seed, linear){
      if(is.matrix(Theta)){
        if(length(linear) == 1){
          X.red <- matrix(X[linear,], nrow =1)
        } else{
          X.red <- X[linear,]
        }
        fitu <- uncollapsePibble_sigmaKnown(coll_samples, X.red, Theta, Gamma, Gamma_comb, Xi, Sigma, upsilon, 
                                            ret_mean, ncores, seed)
        LambdaX <- array(NA, dim = c(nrow(fitu$Lambda), ncol(X), dim(fitu$Lambda)[3]))
        for(j in 1:dim(fitu$Lambda)[3]){
          LambdaX[,,j] <- fitu$Lambda[,,j] %*% X.red
        }
        return(list(Lambda = LambdaX, Lambda.out = fitu$Lambda))
      } else{
        fitu <- uncollapsePibble_sigmaKnown(coll_samples, diag(ncol(X)), Theta(X), Gamma(X), Gamma_comb, Xi, Sigma, upsilon, 
                                            ret_mean, ncores, seed)
        return(list(Lambda = fitu$Lambda, Lambda.out = fitu$Lambda))
      }
    } 
    
    ## setting newdata <- X if newdata is null
    if(is.null(newdata)){
      newdata <- X
    }
    
    for(i in 1:num.comp){
      ## if num.comp == 1 --> return the samples from Lambda above
      if(num.comp == 1){
        ## return the pibble samples.
        if(is.matrix(Theta[[i]])){
          fitu <- uncollapsePibble(collapse_samps$Eta, newdata[linear,], Theta[[i]], Gamma[[i]], Xi, upsilon, 
                                              ret_mean=ret_mean, ncores=ncores, seed=seed)
          Lambda.out[[i]] <- fitu$Lambda
        } else{
          Lambda.out[[i]] <- collapse_samps$Lambda
        }
      } else {
        if(i < num.comp){
          Lambda_mean <- Reduce("+", Lambda)
          Theta_mean <- Reduce("+", Theta_trans[c((i+1):num.comp)])
          unc_samples <- if(is.null(Lambda_mean)){
            sweep(collapse_samps$Lambda, c(1,2), Theta_mean, FUN = "-")
          } else{
            samp_mean <- sweep(Lambda_mean, c(1,2), Theta_mean, FUN = "+")
            collapse_samps$Lambda - samp_mean
          }
          Gamma_comb_red <- Reduce('+', Gamma_trans[c((i+1):num.comp)])
          fitu <- add.uncollapse(unc_samples, newdata, Theta[[i]], Gamma[[i]], Gamma_comb_red, Xi, collapse_samps$Sigma,
                                 upsilon, ret_mean, ncores, seed, linear)
          Lambda[[i]] <- fitu$Lambda
          Lambda.out[[i]] <- fitu$Lambda.out
        } else {
          samp_mean <- Reduce("+", Lambda)
          Lambda[[i]] <- collapse_samps$Lambda - samp_mean
          Lambda.out[[i]] <- Lambda[[i]]
        }
      }
    }

    ## pretty output ##
    out <- list()
    if ("Eta" %in% pars){
      out[["Eta"]] <- collapse_samps$Eta
    }
    if ("Lambda" %in% pars){
      out[["Lambda"]] <- Lambda.out
    }
    if ("Sigma" %in% pars){
      out[["Sigma"]] <- collapse_samps$Sigma
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
    out$linear <- linear
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
    
    ##newdata = NULL automatically handled by pibble
    out <- pibble(Y, X=diag(ncol(X)), upsilon, Theta_train, Gamma_train, Xi, init, pars, newdata = newdata, ...)
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

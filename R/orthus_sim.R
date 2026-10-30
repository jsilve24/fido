#' Simulate simple orthus dataset and priors (for testing)
#' 
#' @param D number of multinomial categories (integer >= 2)
#' @param P number of dimensions of second dataset Z (integer >= 1)
#' @param N number of samples (integer >= 1)
#' @param Q number of covariates (integer >= 1; first one is an intercept)
#' @param use_names should samples, covariates, and categories be named
#' @param true_priors should Xi and upsilon be chosen to have mean at true 
#'   simulated value
#' @return list
#' @export
#' @importFrom stats rnorm rmultinom
#' @examples 
#' sim <- orthus_sim()
orthus_sim <- function(D=10, P=10, N=30, Q=2, use_names=TRUE, true_priors=FALSE){
  dimensions <- list(D = D, P = P, N = N, Q = Q)
  for (name in names(dimensions)) {
    value <- dimensions[[name]]
    minimum <- if (name == "D") 2 else 1
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value != floor(value) || value < minimum) {
      stop(sprintf("%s must be a single integer >= %d", name, minimum))
    }
  }

  # Simulate Data
  Sigma <- diag(sample(1:8, D-1+P, replace=TRUE))
  # A unit off-diagonal is singular when both variances equal one.
  if (D-1+P >= 3 && Sigma[2,2]*Sigma[3,3] > 1) {
    Sigma[2, 3] <- Sigma[3,2] <- -1
  }
  Gamma <- diag(sqrt(rnorm(Q)^2), nrow=Q, ncol=Q)
  Theta <- matrix(0, D-1+P, Q)
  Phi <-  Theta + t(chol(Sigma))%*%matrix(rnorm(Q*(D-1+P)), nrow=D-1+P)%*%chol(Gamma)
  X <- matrix(rnorm(N*(Q-1)), Q-1, N)
  X <- rbind(1, X)
  Eta <- Phi%*%X + t(chol(Sigma))%*%matrix(rnorm(N*(D-1+P)), nrow=D-1+P)
  Z <- Eta[D:(D-1+P), , drop=FALSE]
  Eta <- Eta[1:(D-1), , drop=FALSE]
  Pi <- t(alrInv(t(Eta)))
  Y <- matrix(0, D, N)
  for (i in 1:N) Y[,i] <- rmultinom(1, sample(5000:10000), prob = Pi[,i])
  if (use_names){
    colnames(X) <- colnames(Y) <- paste0("s", 1:N)
    rownames(Y) <- paste0("c", 1:D)
    rownames(X) <- paste0("x", 1:Q)  
  }
  
  # Priors
  if (true_priors){
    upsilon <- D+P+10
    Xi <- Sigma*(upsilon-D-P)
  } else {
    upsilon <- D+P
    Xi <- diag(D-1+P)
  }
  
  # Precompute
  KInv <- solve(Xi)
  AInv <- solve(diag(N)+ t(X)%*%Gamma%*%X)
  
  return(list(Sigma=Sigma, Gamma=Gamma, D=D, N=N, Q=Q, P=P, Theta=Theta, Phi=Phi, 
              X=X, Y=Y,Z=Z, Eta=Eta, upsilon=upsilon, Xi=Xi, KInv=KInv, AInv=AInv))
}

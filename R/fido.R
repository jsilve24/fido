#' fido: Fitting and Analysis of Multinomial Logistic Normal  Models
#' 
#'  Provides methods for fitting and inspection of Bayesian Multinomial 
#'  Logistic Normal Models using MAP estimation and Laplace Approximation.
#'  Key functionality is implemented in C++ for scalability. 
#'  
#' @name fido_package
#' 
#' @useDynLib fido
#' @importFrom Rcpp sourceCpp
"_PACKAGE"

globalVariables(".")


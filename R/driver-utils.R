###Logratio functions

#' Create a default ILR base
#' @param D	the number of parts (e.g., number of columns in untransformed data)
#' @return A matrix
#' @export
create_default_ilr_base <- function(D){
  qr.Q(qr(create_alr_base(D, D)))
}

#' Compute the ALR of a matrix
#'
#' @param x A matrix where the rows are the samples
#' @param d Index of column used as a reference. Defaults to last column
#'
#' @return matrix
#' @export
alr <- function(x, d=NULL){
  x <- vec_to_mat(x)
  if (is.null(d)) d <- ncol(x)
  B <- create_alr_base(ncol(x), d, inv=FALSE)
  glr(x, B)
}

#' Compute the inverse ALR of a matrix
#'
#' @param y An ALR transformed matrix
#' @param d Index of column used as a reference. Defaults to last column
#'
#' @return matrix
#' @export
alrInv <- function(y, d=NULL){
  y <- vec_to_mat(y)
  if (is.null(d)) d <- ncol(y)+1
  B <- create_alr_base(ncol(y)+1, d, inv=TRUE)
  glrInv(y, B)
}

create_alr_base <- function(D,d, inv=FALSE){
  if (d < 1 | d > D) stop("invalid d given D")
  B <- diag(D)
  if (!inv){
    B[d,] <- rep(-1, D)
  } else {
    B[d,] <- rep(0, D)
  }
  B[,-d]
}

glr <- function(x, V){
  log(x) %*% V
}

glrInv <- function(y, V){
  tmp <- exp(y %*% t(V))
  miniclo(tmp)
}

glr_array <- function(x, V, parts, dimname = colnames(V)){
  if (parts == 1){
    dn <- dimnames(x)
    d <- dim(x)
    x <- matrix(x, d[1], prod(d[-1]))
    x <- t(V) %*% log(x)
    d[1] <- ncol(V)
    dim(x) <- d
    if (!is.null(dn)){
      if (!is.null(dimname)){
        dn[[1]] <- dimname
      } else {
        dn[1] <- list(NULL)
      }
      dimnames(x) <- dn
    }
    return(x)
  } else {
    f <- function(x) glr(x, V)
    return(array_apply_1D_function(x, parts, f, dimname))
  }
}

#' Compute the CLR of an array
#'
#' @param x multidimensional array in index
#' @param parts index of dimension of `x` that represents parts
#'
#' @return array
#' @export
clr_array <- function(x, parts){
  n.parts <- dim(x)[parts]
  V <- create_clr_base(n.parts)
  glr_array(x, V, parts)
}

clr <- function(x){
  x <- vec_to_mat(x)
  glr(x, create_clr_base(ncol(x)))
}

create_clr_base <- function(D, inv=FALSE){
  if(!inv){
    M <- matrix(-1, D, D) + D*diag(D)
    return(M/D)
  } else {
    return(diag(D))
  }
}

alrInv_array <- function(y, d=dim(y)[coords]+1, coords){
  B <- create_alr_base(dim(y)[coords]+1, d, inv=TRUE)
  glrInv_array(y, B, coords)
}

glrInv_array <- function(y, V, coords, dimname = rownames(V)){
  if (coords==1){
    dn <- dimnames(y)
    d <- dim(y)
    y <- matrix(y, d[1], prod(d[-1]))
    y <- exp(V %*% y)
    y <- t(miniclo(t(y)))
    d[1] <- nrow(V)
    dim(y) <- d
    if (!is.null(dn)){
      if (!is.null(dimname)) { dn[[1]] <- dimname } else {dn[1] <- list(NULL)}
      dimnames(y) <- dn
    }
    return(y)
  } else {
    f <- function(y) glrInv(y, V)
    return(array_apply_1D_function(y, coords, f, dimname))
  }
}

clrInv_array <- function(y, coords){
  n.coords <- dim(y)[coords]
  V <- diag(n.coords) # Not efficient but reuses code...
  glrInv_array(y, V, coords)
}

ilr_array <- function(x, V=NULL, parts){
  n.parts <- dim(x)[parts]
  if (is.null(V)) V <- create_default_ilr_base(n.parts)
  glr_array(x, V, parts)
}

ilrInv_array <- function(y, V=NULL, coords){
  n.coords <- dim(y)[coords]
  if (is.null(V)) V <- qr.Q(qr(create_alr_base(n.coords+1, n.coords+1)))
  glrInv_array(y, V, coords)
}

#' Compute the ALR of an array
#'
#' @param x multidimensional array in simplex
#' @param d Index of column used as a reference. Defaults to last column
#' @param parts index of dimension of `x` that represents parts
#'
#' @return array
#' @export
alr_array <- function(x, d=dim(x)[parts], parts){
  B <- create_alr_base(dim(x)[parts], d, inv=FALSE)
  glr_array(x, B, parts)
}

#' Compute the ALR of an array
#'
#' @param y multidimensional ALR transformed array
#' @param d Index of column used as a reference. Defaults to last column
#' @param coords index of dimension of `x` that represents coordinates
#'
#' @return array
#' @export
alrInv_array <- function(y, d=dim(y)[coords]+1, coords){
  B <- create_alr_base(dim(y)[coords]+1, d, inv=TRUE)
  glrInv_array(y, B, coords)
}

clrvar2alrvar <- function(Sigma, d2){
  D <- nrow(Sigma)
  V <- create_alr_base(D, d2, inv=FALSE)
  return(t(V)%*%Sigma%*%V)
}

alrvar2clrvar <- function(Sigma, d1){
  D <- nrow(Sigma)+1
  G1 <- create_alr_base(D, d1, inv=TRUE) - 1/D
  return(G1%*%Sigma%*%t(G1))
}

alrvar2alrvar <- function(Sigma, d1, d2){
  S <- alrvar2clrvar(Sigma, d1)
  clrvar2alrvar(S, d2)
}

alrvar2ilrvar <- function(Sigma, d1, V2){
  S <- alrvar2clrvar(Sigma, d1)
  clrvar2ilrvar(S, V2)
}

ilrvar2alrvar <- function(Sigma, V1, d2){
  S <- ilrvar2clrvar(Sigma, V1)
  clrvar2alrvar(S, d2)
}

alrvar2varmat <- function(Sigma, d1){
  S <- alrvar2clrvar(Sigma, d1)
  clrvar2varmat(S)
}

ilrvar2varmat <- function(Sigma, V){
  Sigma <- ilrvar2clrvar(Sigma, V)
  clrvar2varmat(Sigma)
}

clrvar2varmat <- function(Sigma){
  varmat <- matrix(0, nrow(Sigma), ncol(Sigma))
  for (i in 1:dim(Sigma)[1]){
    for (j in 1:dim(Sigma)[2]){
      varmat[i,j] <- Sigma[i,i] + Sigma[j,j] - 2*Sigma[i,j]
    }
  }
}

clrvar2ilrvar <- function(Sigma, V){
  t(V) %*% Sigma %*% V
}

ilrvar2clrvar <- function(Sigma, V){
  V %*% Sigma %*% t(V)
}

ilrvar2ilrvar <- function(Sigma, V1, V2){
  t(V2) %*% V1 %*% Sigma %*% t(V1) %*% V2
}


ilrvar2phi <- function(Sigma, V){
  Sigma.clr <- ilrvar2clrvar(Sigma,V)
  phi <- Sigma.clr
  phi[] <- 0
  for(i in 1:dim(Sigma.clr)[1]){
    for (j in 1:dim(Sigma.clr)[2]){
      phi[i,j] <- Sigma.clr[i,i] + Sigma.clr[j,j] - 2*Sigma.clr[i,j]
    }
  }
  return(phi)
}


##Helper functions

#' Gather Multidimensional Array to Tidy Tibble
#'
#' @param a multidimensional array
#' @param value unquoted name of column with values (defaults to "var")
#' @param ... unquoted dimension names (defaults to "dim_1", "dim_2", etc...)
#' @param .id if specified, name for column created with name of a captured
#'
#' @return data.frame
#' @seealso spread_array
#' @export
#' @import dplyr purrr tidyr
#' @importFrom rlang quos enquo quo_name sym syms
#'
#' @examples
#' a <- array(1:100, dim =c(10, 5, 2))
#' gather_array(a, sequence, A, B, C)
gather_array <- function(a, value, ..., .id=NULL){
  qs <- rlang::quos(...)
  if (missing(value)) {
    evalue <- rlang::sym("var")}
  else {
    evalue <- rlang::enquo(value)
  }
  len <- length(qs)
  d <- dim(a)
  
  # Default Values
  if (len > 0) {
    dimnames <- purrr::map(qs, rlang::quo_name) %>%
      purrr::as_vector()
  } else {
    dimnames <- paste0("dim_", 1:length(d))
  }
  
  l <- list()
  for (i in 1:length(d)){
    l[[i]] <- 1:d[i]
  }
  names(l) <- dimnames
  tidy <- expand.grid(l)
  tidy[[rlang::quo_name(evalue)]] <- a[as.matrix(tidy)]
  if (!is.null(.id)) tidy[[.id]] <- rlang::expr_name(a)
  return(tidy)
}

#' Shortcut for summarize variable with quantiles and mean
#'
#' @param data tidy data frame
#' @param var variable name (unquoted) to be summarised
#' @param ... other expressions to pass to summarise
#'
#' @return data.frame
#' @export
#' @details Notation: \code{pX} refers to the \code{X}\% quantile
#' @import dplyr
#' @importFrom stats quantile
#' @importFrom rlang quos quo UQ
#' @examples
#' d <- data.frame("a"=sample(1:10, 50, TRUE),
#'                 "b"=rnorm(50))
#'
#' # Summarize posterior for b over grouping of a and also calcuate
#' # minmum of b (in addition to normal statistics returned)
#' d <- dplyr::group_by(d, a)
#' summarise_posterior(d, b, mean.b = mean(b), min=min(b))
summarise_posterior <- function(data, var, ...){
  qvar <- enquo(var)
  qs <- quos(...)
  
  
  data %>%
    summarise(p2.5 = quantile(!!qvar, prob=0.025),
              p25 = quantile(!!qvar, prob=0.25),
              p50 = quantile(!!qvar, prob=0.5),
              mean = mean(!!qvar),
              p75 = quantile(!!qvar, prob=0.75),
              p97.5 = quantile(!!qvar, prob=0.975),
              !!!qs)
}

vec_to_mat <- function(x){
  if (is.vector(x)) {
    n <- names(x)
    x <- matrix(x, nrow = 1)
    colnames(x) <- n
  }
  x
}

#' @importFrom rlang :=
array_apply_1D_function <- function(a, dimno, f, dimname=NULL){
  
  d <- dim(a)
  ndim <- length(d)
  sdim <- sym(paste0("dim_", dimno))
  sdim_other <- syms(paste0("dim_", (1:ndim)[(1:ndim) != dimno]))
  
  # Store Dimnames
  dn <- dimnames(a)
  
  # Actual Computation
  var <- "var"
  ga <- a %>%
    gather_array(!!var) %>%
    spread(!!sdim, var)
  
  indicies <- ga %>%
    select(!!!sdim_other)
  
  b <- ga %>%
    select(-contains("dim")) %>%
    as.matrix() %>%
    f %>%
    `colnames<-`(., 1:ncol(.)) %>%
    as.data.frame() %>%
    bind_cols(indicies, .) %>%
    gather(!!sdim, var, -contains("dim")) %>%
    mutate(!!quo_name(sdim)  := as.integer(!!sdim)) %>%
    spread_array(var, !!!syms(paste0("dim_", 1:ndim)))
  
  # Update Dimnames
  if (!is.null(dn) || !is.null(dimname)){
    if (!is.null(dn)){
      dn[dimno] <- list(NULL)
    } else {
      dn <- list()
      for (i in 1:length(d)) dn[[i]] <- NULL
    }
    if (!is.null(dimname)){
      dn[[dimno]] <- dimname
    }
    dimnames(b) <- dn
  } else {
    names(dim(b)) <- NULL
  }
  return(b)
}

add_array_dim <- function(a, d){
  dd <- dim(a)
  if (d > length(dd)+1) stop("d must be <= length(dim(a))+1")
  ad <- rep(NA, length(dd)+1)
  passed=FALSE
  for (i in 1:(length(dd)+1)) {
    if (i==d) { ad[i] <- 1; passed <- TRUE }
    else if (passed) { ad[i] <- dd[i-1] }
    else {ad[i] <- dd[i]}
  }
  array(a, dim=ad)
}

spread_array <- function(data, value, ...){
  evalue <- rlang::enquo(value)
  qs <- rlang::quos(...)
  l <- length(qs)
  
  # Default Values
  if (l == 0) {
    cn <- colnames(data)
    cn <- cn[grepl("dim_", cn)]
    # Validation of Defaults
    consecutive1 <- strsplit(cn, "_") %>%
      purrr::map(~.x[2]) %>%
      purrr::map(as.integer) %>%
      as_vector()
    if (!setequal(consecutive1, 1:length(consecutive1))) {
      stop("default dimnetion names must have consecutive integer suffixes")
    }
    ##
    cn <- cn[match(consecutive1, 1:length(consecutive1))]
    qs <- rlang::syms(cn)
  }
  if (missing(value)) evalue <- rlang::sym("var")
  
  tidy_dim <- data %>%
    select(!!!qs)
  unique_dim <- tidy_dim %>%
    as.list() %>%
    purrr::map(unique)
  length_dim <- unique_dim %>%
    purrr::map(length)
  
  # Input validation - Must be sequential integers
  class_dim <- data %>%
    select(!!!qs) %>%
    sapply(class)
  if (!all(class_dim=="integer")) stop("Dimension indexes must be integers")
  
  consecutive2 <- unique_dim %>%
    map2(map(length_dim, ~1:.x), setequal) %>%
    as_vector() %>%
    all()
  if (!consecutive2) stop("Dimension indexes must be consecutive")
  #####
  
  a <- array(NA, dim = length_dim)
  a[as.matrix(tidy_dim)] <- pull(data, rlang::quo_name(evalue))
  
  return(a)
}

##Closure functions

#' Closure operator
#'
#' @param x vector or matrix (rows are samples, parts are columns) of data in simplex
#'
#' @return x with row entries divided by sum of row (converts vectors to row matricies)
#' @export
#'
#' @examples
#' x <- matrix(runif(30), 10, 3)
#' x <- miniclo(x)
miniclo <- function(x){
  if (is.vector(x)) {
    n <- names(x)
    x <- matrix(x, nrow = 1)
    colnames(x) <- n
  }
  (x/rowSums(x))
}

#' Closure Operation applied to array on margin
#'
#' Array version of \code{\link{miniclo}}.
#'
#' @param x multidimensional array
#' @param parts index of dimension of \code{x} that represents parts (e.g., compositional variables)
#'
#' @return array
#' @export
#'
#' @examples
#' x <- array(1:100, dim=c(10, 5, 2))
#' miniclo_array(x, parts=2)
miniclo_array <- function(x, parts){
  array_apply_1D_function(x, parts, miniclo)
}
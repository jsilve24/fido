#' Transform Fit fido Parameters to other representations
#' 
#' These are a collection of convenience functions for transforming
#' fido fit objects to a number of different representations including
#' ILR bases, CLR coordinates, ALR coordinates, and proportions. 
#' 
#' For orthus, transforms only appleid to log-ratio parameters
#' 
#' @param m object of class pibblefit or orthusfit (e.g., output of \code{\link{pibble}}
#'   or \code{\link{orthus}})
#' @param d (integer) multinomial category to take as new alr reference
#' @param V (matrix) contrast matrix for ILR basis to transform into to (defaults to 
#'   \code{create_default_ilr_base(D)})
#'
#' @details Note: that there is a degeneracy of representations for a covariance 
#' matrix represented in terms of proportions. As such the function 
#' \code{to_proportions} does not attempt to transform parameters Sigma
#' or prior Xi and instead just removes them from the pibblefit object returned. 
#' 
#' @return object
#' @name fido_transforms
NULL


#' @rdname fido_transforms
#' @export
to_proportions <- function(m){
  UseMethod("to_proportions",m)
}

#' @rdname fido_transforms
#' @export
to_alr <- function(m, d){
  UseMethod("to_alr",m)
}

#' @rdname fido_transforms
#' @export
to_ilr <- function(m, V=NULL){
  UseMethod("to_ilr", m)
}

#' @rdname fido_transforms
#' @export
to_clr <- function(m){
  UseMethod("to_clr",m)
}

#' @rdname fido_transforms
#' @export
to_proportions.pibblefit <- function(m){
  if (m$coord_system == "alr"){
    if (!is.null(m$Eta)) m$Eta <- alrInv_array(m$Eta, m$alr_base, 1)
    if(!is.null(m$Lambda)){
      if(is.list(m$Lambda)){
        for(i in 1:length(m$Lambda)) m$Lambda[[i]] <-  alrInv_array(m$Lambda[[i]], m$alr_base, 1)
      }
      else {
        m$Lambda <-  alrInv_array(m$Lambda, m$alr_base, 1)
      }
    }
    if (!is.null(m$Sigma)){
      if (m$alr_base != m$D){
        for (i in 1:m$iter){
          m$Sigma[,,i] <- alrvar2alrvar(m$Sigma[,,i], m$alr_base, m$D)
        }
      }
      m$Sigma_default <- m$Sigma
      m$Sigma <- NULL
    }
    # Transform Priors as well 
    if (!is.null(m$Xi)){
      if (m$alr_base != m$D){
        m$Xi <- alrvar2alrvar(m$Xi, m$alr_base, m$D)
      }
      m$Xi_default <- m$Xi
      m$Xi <- NULL
    }
    if (!is.null(m$Theta)){
      if (!inherits(m, "bassetfit")) m$Theta <- alrInv_array(m$Theta, m$alr_base, 1)
    }
    if (!is.null(m$init)) m$init <- alrInv_array(m$init, m$alr_base, 1)
  }
  if (m$coord_system == "ilr"){
    if (!is.null(m$Eta)) m$Eta <- ilrInv_array(m$Eta, m$ilr_base, 1)
    if(!is.null(m$Lambda)){
      if(is.list(m$Lambda)){
        for(i in 1:length(m$Lambda)) m$Lambda[[i]] <-  ilrInv_array(m$Lambda[[i]], m$ilr_base, 1)
      }
      else {
        m$Lambda <-  ilrInv_array(m$Lambda, m$ilr_base, 1)
      }
    }
    if (!is.null(m$Sigma)){
      for (i in 1:m$iter){
        m$Sigma[,,i] <- ilrvar2alrvar(m$Sigma[,,i], m$ilr_base, m$D)
      }
      m$Sigma_default <- m$Sigma
      m$Sigma <- NULL
    }
    
    # Transform priors as well 
    if (!is.null(m$Xi)){
      m$Xi <- ilrvar2alrvar(m$Xi, m$ilr_base, m$D)
      m$Xi_default <- m$Xi
      m$Xi <- NULL  
    }
    if (!is.null(m$Theta)) {
      if (!inherits(m, "bassetfit")) m$Theta <- ilrInv_array(m$Theta, m$ilr_base, 1)  
    }
    if (!is.null(m$init)) m$init <- ilrInv_array(m$init, m$ilr_base, 1)
  }
  if (m$coord_system == "clr"){
    if (!is.null(m$Eta)) m$Eta <- clrInv_array(m$Eta, 1)
    if(!is.null(m$Lambda)){
      if(is.list(m$Lambda)){
        for(i in 1:length(m$Lambda)) m$Lambda[[i]] <- clrInv_array(m$Lambda[[i]], 1)
      }
      else {
        m$Lambda <- clrInv_array(m$Lambda, 1)
      }
    }
    if (!is.null(m$Sigma)){
      Sigma_default <- array(0, dim=c(m$D-1, m$D-1, m$iter))
      for (i in 1:m$iter){
        Sigma_default[,,i] <- clrvar2alrvar(m$Sigma[,,i], m$D)
      }
      m$Sigma <- NULL
      m$Sigma_default <- Sigma_default
    }
    # Transform priors as well
    if (!is.null(m$Xi)){
      m$Xi_default <- clrvar2alrvar(m$Xi, m$D)
      m$Xi <- NULL      
    }
    if (!is.null(m$Theta)){
      if (!inherits(m, "bassetfit")) m$Theta <- clrInv_array(m$Theta, 1)  
    }
    if (!is.null(m$init)) m$init <- clrInv_array(m$init, 1)
  }
  if (m$coord_system=="proportions"){
    return(m)
  }
  m$summary <- NULL
  m$coord_system <- "proportions"
  m$ilr_base <- NULL
  m$alr_base <- NULL
  return(m)
}


#' @rdname fido_transforms
#' @export
to_proportions.orthusfit <- function(m){
  if (m$coord_system == "alr"){
    if (!is.null(m$Eta)) m$Eta <- alrInv_array(m$Eta, m$alr_base, 1) 
    if (!is.null(m$Lambda)) m$Lambda <- oalrInv(m$Lambda, m$D-1, m$alr_base)
    
    if (!is.null(m$Sigma)){
      if (m$alr_base != m$D) m$Sigma <- oalrvar2alrvar(m$Sigma, m$D-1, m$alr_base, m$D) 
      m$Sigma_default <- m$Sigma
      m$Sigma <- NULL
    }
    # Transform Priors as well  
    if (!is.null(m$Xi)){
      if (m$alr_base != m$D) m$Xi <- oalrvar2alrvar(m$Xi, m$D-1, m$alr_base, m$D)
      m$Xi_default <- m$Xi
      m$Xi <- NULL
    }
    if (!is.null(m$Theta)) m$Theta <- oalrInv(m$Theta, m$D-1, m$alr_base)
    if (!is.null(m$init)) m$init <- alrInv_array(m$init, m$alr_base, 1) 
  }
  if (m$coord_system == "ilr"){
    if (!is.null(m$Eta)) m$Eta <- ilrInv_array(m$Eta, m$ilr_base, 1)
    if (!is.null(m$Lambda)) m$Lambda <- oilrInv(m$Lambda, m$D-1, m$ilr_base)
    if (!is.null(m$Sigma)) {
      m$Sigma_default <- oilrvar2alrvar(m$Sigma, m$D-1, m$ilr_base, m$D) 
      m$Sigma <- NULL
    }
    
    # Transform priors as well 
    if (!is.null(m$Xi)){
      m$Xi_default <- oilrvar2alrvar(add_array_dim(m$Xi,3), m$D-1, m$ilr_base, m$D)[,,1]
      m$Xi <- NULL  
    }
    if (!is.null(m$Theta)) m$Theta <- oilrInv(m$Theta, m$D-1, m$ilr_base)
    if (!is.null(m$init)) m$init <- ilrInv_array(m$init, m$ilr_base, 1)
  }
  if (m$coord_system == "clr"){
    if (!is.null(m$Eta)) m$Eta <- clrInv_array(m$Eta, 1)
    if (!is.null(m$Lambda)) m$Lambda <- oclrInv(m$Lambda, m$D)
    if (!is.null(m$Sigma)){
      m$Sigma_default <- oclrvar2alrvar(m$Sigma, m$D, m$D)
      m$Sigma <- NULL
    }
    # Transform priors as well
    if (!is.null(m$Xi)){
      m$Xi_default <- oclrvar2alrvar(m$Xi, m$D, m$D)
      m$Xi <- NULL      
    }
    if (!is.null(m$Theta)) m$Theta <- oclrInv(m$Theta, m$D)  
    if (!is.null(m$init)) m$init <- clrInv_array(m$init,1)
  }
  if (m$coord_system=="proportions"){
    return(m)
  }
  m$summary <- NULL
  m$coord_system <- "proportions"
  m$ilr_base <- NULL
  m$alr_base <- NULL
  return(m)
}

#' @rdname fido_transforms
#' @export
to_alr.pibblefit <- function(m, d){
  if (m$coord_system=="alr"){
    if (m$alr_base == d) return(m)
  }
  m <- to_proportions(m)
  
  if (!is.null(m$Eta)) m$Eta <- alr_array(m$Eta, d, 1)
  
  if(!is.null(m$Lambda)){
    if(is.list(m$Lambda)){
      for(i in 1:length(m$Lambda)) m$Lambda[[i]] <- alr_array(m$Lambda[[i]], d, 1)
    }
    else {
      m$Lambda <- alr_array(m$Lambda, d, 1)
    }
  }
  if (!is.null(m$Sigma)){
    m$Sigma <- array(0, dim=dim(m$Sigma_default))
    for (i in 1:m$iter){
      m$Sigma[,,i] <- alrvar2alrvar(m$Sigma_default[,,i], m$D, d)
    }
    m$Sigma_default <- NULL
  }
  # Transform priors as well 
  if (!is.null(m$Xi)){
    m$Xi <- alrvar2alrvar(m$Xi_default, m$D, d)
    m$Xi_default <- NULL  
  }
  if (!is.null(m$Theta)){
    if (!inherits(m, "bassetfit")) m$Theta <- alr_array(m$Theta, d, 1)  
  }
  if (!is.null(m$init)) m$init <- alr_array(m$init, d, 1)
  
  m$summary <- NULL
  m$coord_system <- "alr"
  m$alr_base <- d
  return(m)
}

#' @rdname fido_transforms
#' @export
to_alr.orthusfit <- function(m, d){
  if (m$coord_system=="alr"){
    if (m$alr_base == d) return(m)
  }
  m <- to_proportions.orthusfit(m)
  
  if (!is.null(m$Eta)) m$Eta <- alr_array(m$Eta, d, 1)
  if (!is.null(m$Lambda)) m$Lambda <- oalr(m$Lambda, m$D, d)
  if (!is.null(m$Sigma)){
    m$Sigma <- oalrvar2alrvar(m$Sigma_default, m$D-1, m$D, d)
    m$Sigma_default <- NULL
  }
  # Transform priors as well 
  if (!is.null(m$Xi)){
    m$Xi <- oalrvar2alrvar(m$Xi_default, m$D-1, m$D, d)
    m$Xi_default <- NULL  
  }
  if (!is.null(m$Theta)) m$Theta <- oalr(m$Theta, m$D, d)  
  if (!is.null(m$init)) m$init <- alr_array(m$init, d, 1)
  
  m$summary <- NULL
  m$coord_system <- "alr"
  m$alr_base <- d
  return(m)
}

#' @rdname fido_transforms
#' @export
to_ilr.pibblefit <- function(m, V=NULL){
  if (m$coord_system=="ilr"){
    if (all.equal(m$ilr_base, V)) return(m)
  }
  if (is.null(V)) V <- create_default_ilr_base(m$D)
  m <- to_proportions(m)
  
  if (!is.null(m$Eta)) m$Eta <- ilr_array(m$Eta, V, 1)
  if(!is.null(m$Lambda)){
    if(is.list(m$Lambda)){
      for(i in 1:length(m$Lambda)) m$Lambda[[i]] <-  ilr_array(m$Lambda[[i]], V, 1)
    }
    else {
      m$Lambda <-  ilr_array(m$Lambda, V, 1)
    }
  }
  if (!is.null(m$Sigma)){
    m$Sigma <- array(0, dim=dim(m$Sigma_default))
    for (i in 1:m$iter){
      m$Sigma[,,i] <- alrvar2ilrvar(m$Sigma_default[,,i], m$D, V)
    }
    m$Sigma_default <- NULL
  }
  # Transform priors as well 
  if (!is.null(m$Xi)){
    m$Xi <- alrvar2ilrvar(m$Xi_default, m$D, V)
    m$Xi_default <- NULL  
  }
  if (!is.null(m$Theta)){
    if (!inherits(m, "bassetfit")) m$Theta <- ilr_array(m$Theta, V, 1)  
  }
  if (!is.null(m$init)) m$init <- ilr_array(m$init, V, 1)
  
  colnames(V) <- paste0("ilr_", rep(1:ncol(V)))
    
  m$summary <- NULL
  m$coord_system <- "ilr"
  m$ilr_base <- V
  return(m)
}

#' @rdname fido_transforms
#' @export
to_ilr.orthusfit <- function(m, V=NULL){
  if (m$coord_system=="ilr"){
    if (all.equal(m$ilr_base, V)) return(m)
  }
  if (is.null(V)) V <- create_default_ilr_base(m$D)
  m <- to_proportions(m)
  
  if (!is.null(m$Eta)) m$Eta <- ilr_array(m$Eta, V, 1)
  if (!is.null(m$Lambda)) m$Lambda <- oilr(m$Lambda, m$D, V)
  if (!is.null(m$Sigma)){
    m$Sigma <- oalrvar2ilrvar(m$Sigma_default, m$D-1, m$D, V)
    m$Sigma_default <- NULL
  }
  # Transform priors as well 
  if (!is.null(m$Xi)){
    m$Xi <- oalrvar2ilrvar(m$Xi_default, m$D-1, m$D, V)
    m$Xi_default <- NULL  
  }
  if (!is.null(m$Theta)) m$Theta <- oilr(m$Theta, m$D, V)  
  if (!is.null(m$init)) m$init <- ilr_array(m$init, V, 1)
  
  colnames(V) <- paste0("ilr_", rep(1:ncol(V)))
  
  m$summary <- NULL
  m$coord_system <- "ilr"
  m$ilr_base <- V
  return(m)
}

#' @rdname fido_transforms
#' @export
to_clr.pibblefit <- function(m){
  if (m$coord_system=="clr") return(m)
  m <- to_proportions(m)

  if (!is.null(m$Eta)) m$Eta <- clr_array(m$Eta, 1)
  if(!is.null(m$Lambda)){
    if(is.list(m$Lambda)){
      for(i in 1:length(m$Lambda)) m$Lambda[[i]] <- clr_array(m$Lambda[[i]], 1)
      }
    else {
      m$Lambda <- clr_array(m$Lambda, 1)
    }
  }
  if (!is.null(m$Sigma)){
    m$Sigma <- array(0, dim=c(m$D, m$D, m$iter))
    for (i in 1:m$iter){
      m$Sigma[,,i] <- alrvar2clrvar(m$Sigma_default[,,i], m$D)
    }
    m$Sigma_default <- NULL
  }
  # Transform priors as well 
  if (!is.null(m$Xi)){
    m$Xi <- alrvar2clrvar(m$Xi_default, m$D)
    m$Xi_default <- NULL  
  }
  if (!is.null(m$Theta)){
    if (!inherits(m, "bassetfit")) m$Theta <- clr_array(m$Theta, 1)  
  }
  if (!is.null(m$init)) m$init <- clr_array(m$init, 1)
  
  m$summary <- NULL
  m$coord_system <- "clr"
  return(m)
}


#' @rdname fido_transforms
#' @export
to_clr.orthusfit <- function(m){
  if (m$coord_system=="clr") return(m)
  m <- to_proportions(m)
  
  if (!is.null(m$Eta)) m$Eta <- clr_array(m$Eta, 1)
  if (!is.null(m$Lambda)) m$Lambda <- oclr(m$Lambda, m$D-1)
  if (!is.null(m$Sigma)){
    m$Sigma <- oalrvar2clrvar(m$Sigma_default, m$D-1, m$D)
    m$Sigma_default <- NULL
  }
  # Transform priors as well 
  if (!is.null(m$Xi)){
    m$Xi <- oalrvar2clrvar(m$Xi_default, m$D-1, m$D)
    m$Xi_default <- NULL  
  }
  if (!is.null(m$Theta)) m$Theta <- oclr(m$Theta, m$D-1)  
  if (!is.null(m$init)) m$init <- clr_array(m$init, 1)
  
  m$summary <- NULL
  m$coord_system <- "clr"
  return(m)
}


#' Holds information on coordinates system to later be reapplied
#' 
#' \code{store_coord} stores coordinate information for pibblefit object
#' and can be reapplied with function \code{reapply_coord}. Some coordinate
#' systems are not useful for computation and this makes it simple keep 
#' returned object from computations in the same coordinate system as the input. 
#' (Likely most useful inside of a package)
#' 
#' 
#' @param m object of class pibblefit
#' @param l object returned by function \code{store_coord}
#' @name store_coord
#' @return \code{store_coord} list with important information to identify c
#'  coordinate system of pibblefit object. \code{reapply_coord} pibblefit object
#'  in coordinate system previously stored. 
NULL


#' @rdname store_coord
#' @export
store_coord <- function(m){
  l <- list()
  l$coord_system <- m$coord_system
  l$alr_base <- m$alr_base
  l$ilr_base <- m$ilr_base
  return(l)
}

#' @rdname store_coord
#' @export
reapply_coord <- function(m, l){
  if (l$coord_system == "proportions") return(to_proportions(m))
  if (l$coord_system == "clr") return(to_clr(m))
  if (l$coord_system == "alr") return(to_alr(m, l$alr_base))
  if (l$coord_system == "ilr") return(to_ilr(m, l$ilr_base))
  stop("not a recognized coordinate system")
}
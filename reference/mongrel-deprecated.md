# mongrel

This function is deprecated, please use `pibble` instead.

## Usage

``` r
mongrel(
  Y = NULL,
  X = NULL,
  upsilon = NULL,
  Theta = NULL,
  Gamma = NULL,
  Xi = NULL,
  init = NULL,
  pars = c("Eta", "Lambda", "Sigma"),
  ...
)
```

## Arguments

- Y:

  D x N matrix of counts (if NULL uses priors only)

- X:

  Q x N matrix of covariates (design matrix) (if NULL uses priors only,
  must be present to sample Eta)

- upsilon:

  dof for inverse wishart prior (numeric must be \> D) (default: D+3)

- Theta:

  (D-1) x Q matrix of prior mean for regression parameters (default:
  matrix(0, D-1, Q))

- Gamma:

  QxQ prior covariance matrix (default: diag(Q))

- Xi:

  (D-1)x(D-1) prior covariance matrix (default: ALR transform of
  diag(1)\*(upsilon-D)/2 - this is essentially iid on "base scale" using
  Aitchison terminology)

- init:

  (D-1) x N initialization for Eta for optimization

- pars:

  character vector of posterior parameters to return

- ...:

  arguments passed to
  [`optimPibbleCollapsed`](https://jsilve24.github.io/fido/reference/optimPibbleCollapsed.md)
  and
  [`uncollapsePibble`](https://jsilve24.github.io/fido/reference/uncollapsePibble.md)

## Value

An object of class pibblefit

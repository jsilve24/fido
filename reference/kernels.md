# Multivariate RBF Kernel

Designed to be partially specified. (see examples)

## Usage

``` r
SE(X, sigma = 1, rho = median(as.matrix(dist(t(X)))), jitter = 1e-10)

LINEAR(X, sigma = 1, c = rep(0, nrow(X)))
```

## Arguments

- X:

  covariate (dimension Q x N; i.e., covariates x samples)

- sigma:

  scalar parameter

- rho:

  scalar bandwidth parameter

- jitter:

  small scalar to add to off-diagonal of gram matrix (for numerical
  underflow issues)

- c:

  vector parameter defining intercept for linear kernel

## Value

Gram Matrix (N x N) (e.g., the Kernel evaluated at each pair of points)

## Details

Gram matrix G is given by

SE (squared exponential): \$\$G = \sigma^2 \*
exp(-\[(X-c)'(X-c)\]/(s\*\rho^2))\$\$

LINEAR: \$\$G = \sigma^2\*(X-c)'(X-c)\$\$

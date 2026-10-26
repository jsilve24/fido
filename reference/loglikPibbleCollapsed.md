# Calculations for the Collapsed Pibble Model

Functions providing access to the Log Likelihood, Gradient, and Hessian
of the collapsed pibble model. Note: These are convenience functions but
are not as optimized as direct coding of the PibbleCollapsed C++ class
due to a lack of Memoization. By contrast function optimPibbleCollapsed
is much more optimized and massively cuts down on repeated calculations.
A more efficient Rcpp module based implementation of these functions may
following if the future. For model details see
[`optimPibbleCollapsed`](https://jsilve24.github.io/fido/reference/optimPibbleCollapsed.md)
documentation. The log likelihood omits terms constant in eta (including
the matrix-t and multinomial normalizers); it is not a normalized
density.

## Usage

``` r
loglikPibbleCollapsed(Y, upsilon, ThetaX, KInv, AInv, eta, sylv = FALSE)

gradPibbleCollapsed(Y, upsilon, ThetaX, KInv, AInv, eta, sylv = FALSE)

hessPibbleCollapsed(Y, upsilon, ThetaX, KInv, AInv, eta, sylv = FALSE)
```

## Arguments

- Y:

  D x N matrix of counts

- upsilon:

  matrix-t degrees of freedom (must be \> 0). This is *not* the same
  scale as the `upsilon` argument of
  [`pibble`](https://jsilve24.github.io/fido/reference/pibble_fit.md) or
  [`uncollapsePibble`](https://jsilve24.github.io/fido/reference/uncollapsePibble.md):
  for an inverse-Wishart prior with degrees of freedom `upsilon`, pass
  `upsilon - D + 2` here.

- ThetaX:

  D-1 x N matrix formed by Theta\*X (Theta is Prior mean for regression
  coefficients)

- KInv:

  Inverse of K for LTP (for Pibble defined as KInv = solve(Xi))

- AInv:

  Inverse of A for LTP (for Pibble defined as AInv = solve(diag(N)+
  X'GammaX) )

- eta:

  matrix (D-1)xN of parameter values at which to calculate quantities

- sylv:

  (default:false) if true and if N \< D-1 will use sylvester determinant
  identity to speed computation

## Value

see below

- loglikPibbleCollapsed - double

- gradPibbleCollapsed - vector

- hessPibbleCollapsed- matrix

## Examples

``` r
D <- 10
Q <- 2
N <- 30

# Simulate Data
Sigma <- diag(sample(1:8, D-1, replace=TRUE))
Sigma[2, 3] <- Sigma[3,2] <- -1
Gamma <- diag(sqrt(rnorm(Q)^2))
Theta <- matrix(0, D-1, Q)
Phi <-  Theta + t(chol(Sigma))%*%matrix(rnorm(Q*(D-1)), nrow=D-1)%*%chol(Gamma)
X <- matrix(rnorm(N*(Q-1)), Q-1, N)
X <- rbind(1, X)
Eta <- Phi%*%X + t(chol(Sigma))%*%matrix(rnorm(N*(D-1)), nrow=D-1)
Pi <- t(alrInv(t(Eta)))
Y <- matrix(0, D, N)
for (i in 1:N) Y[,i] <- rmultinom(1, sample(5000:10000), prob = Pi[,i])

# Priors
upsilon <- D+10
Xi <- Sigma*(upsilon-D)

# Precompute
KInv <- solve(Xi)
AInv <- solve(diag(N)+ t(X)%*%Gamma%*%X)
ThetaX <- Theta%*%X


loglikPibbleCollapsed(Y, upsilon - D + 2, ThetaX, KInv, AInv, Eta)
#> [1] -195051.6
gradPibbleCollapsed(Y, upsilon - D + 2, ThetaX, KInv, AInv, Eta)[1:5]
#> [1] -4.6504749 -0.2351229  1.2316352  0.2345338  2.8531653
hessPibbleCollapsed(Y, upsilon - D + 2, ThetaX, KInv, AInv, Eta)[1:5,1:5]
#>               [,1]        [,2]        [,3]        [,4]          [,5]
#> [1,] -395.52246082  0.39157593 -0.01261465 -0.05799243  393.92071965
#> [2,]    0.39157593 -8.25804838 -0.26758007 -0.09229975    6.34270080
#> [3,]   -0.01261465 -0.26758007 -1.17258369  0.08725493    0.80984797
#> [4,]   -0.05799243 -0.09229975  0.08725493 -0.42836088    0.04745915
#> [5,]  393.92071965  6.34270080  0.80984797  0.04745915 -413.33444034
```

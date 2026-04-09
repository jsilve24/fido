# Calculations for the Collapsed Pibble Model

Functions providing access to the Log Likelihood, Gradient, and Hessian
of the collapsed pibble model. Note: These are convenience functions but
are not as optimized as direct coding of the PibbleCollapsed C++ class
due to a lack of Memoization. By contrast function optimPibbleCollapsed
is much more optimized and massively cuts down on repeated calculations.
A more efficient Rcpp module based implementation of these functions may
following if the future. For model details see
[`optimPibbleCollapsed`](https://jsilve24.github.io/fido/reference/optimPibbleCollapsed.md)
documentation

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

  (must be \> D)

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


loglikPibbleCollapsed(Y, upsilon, ThetaX, KInv, AInv, Eta)
#> [1] -195095.9
gradPibbleCollapsed(Y, upsilon, ThetaX, KInv, AInv, Eta)[1:5]
#> [1] -4.7284477 -0.1509676  1.2405685  0.2829217  2.6895085
hessPibbleCollapsed(Y, upsilon, ThetaX, KInv, AInv, Eta)[1:5,1:5]
#>               [,1]       [,2]        [,3]        [,4]          [,5]
#> [1,] -395.55832129  0.3806120 -0.02360702 -0.06801042  393.92582730
#> [2,]    0.38061199 -8.4987346 -0.31053592 -0.10707949    6.35296890
#> [3,]   -0.02360702 -0.3105359 -1.22834285  0.10121428    0.81696424
#> [4,]   -0.06801042 -0.1070795  0.10121428 -0.48603619    0.04496538
#> [5,]  393.92582730  6.3529689  0.81696424  0.04496538 -413.36822195
```

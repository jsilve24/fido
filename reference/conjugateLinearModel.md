# Solve Bayesian Multivariate Conjugate Linear Model

See details for model. Notation: `N` is number of samples, `D` is the
dimension of the response, `Q` is number of covariates.

## Usage

``` r
conjugateLinearModel(Y, X, Theta, Gamma, Xi, upsilon, n_samples = 2000L)
```

## Arguments

- Y:

  matrix of dimension D x N

- X:

  matrix of covariates of dimension Q x N

- Theta:

  matrix of prior mean of dimension D x Q

- Gamma:

  covariance matrix of dimension Q x Q

- Xi:

  covariance matrix of dimension D x D

- upsilon:

  scalar (must be \> D) degrees of freedom for InvWishart prior

- n_samples:

  number of samples to draw (default: 2000)

## Value

List with components

1.  Lambda Array of dimension D x Q x n_samples (posterior samples)

2.  Sigma Array of dimension D x D x n_samples (posterior samples)

## Details

\$\$Y \sim MN\_{D \times N}(\Lambda \mathbf{X}, \Sigma, I_N)\$\$
\$\$\Lambda \sim MN\_{D \times Q}(\Theta, \Sigma, \Gamma)\$\$ \$\$\Sigma
\sim InvWish(\upsilon, \Xi)\$\$ This function provides a means of
sampling from the posterior distribution of `Lambda` and `Sigma`.

## Examples

``` r
sim <- pibble_sim()
eta.hat <- t(alr(t(sim$Y+0.65)))
fit <- conjugateLinearModel(eta.hat, sim$X, sim$Theta, sim$Gamma, 
                            sim$Xi, sim$upsilon, n_samples=2000)
```

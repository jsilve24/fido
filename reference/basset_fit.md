# Interface to fit basset models

Basset (A Lazy Learner) - non-linear regression models in fido

## Usage

``` r
basset(
  Y = NULL,
  X,
  upsilon = NULL,
  Theta = NULL,
  Gamma = NULL,
  Xi = NULL,
  linear = NULL,
  init = NULL,
  pars = c("Eta", "Lambda", "Sigma"),
  newdata = NULL,
  ...
)

# S3 method for class 'bassetfit'
refit(m, pars = c("Eta", "Lambda", "Sigma"), ...)
```

## Arguments

- Y:

  D x N matrix of counts (if NULL uses priors only)

- X:

  Q x N matrix of covariates (cannot be NULL)

- upsilon:

  dof for inverse wishart prior (numeric must be \> D) (default: D+3)

- Theta:

  A function from dimensions dim(X) -\> (D-1)xN (prior mean of gaussian
  process). For an additive GP model, can be a list of functions from
  dimensions dim(X) -\> (D-1)xN + a (optional) matrix of size (D-1)xQ
  for the prior of a linear component if desired.

- Gamma:

  A function from dimension dim(X) -\> NxN (kernel matrix of gaussian
  process). For an additive GP model, can be a list of functions from
  dimension dim(X) -\> NxN + a QxQ prior covariance matrix if a linear
  component is specified. It is assumed that the order matches the order
  of Theta.

- Xi:

  (D-1)x(D-1) prior covariance matrix (default: ALR transform of
  diag(1)\*(upsilon-D)/2 - this is essentially iid on "base scale" using
  Aitchison terminology)

- linear:

  A vector denoting which rows of X should be used if a linear component
  is specified. Default is all rows.

- init:

  (D-1) x Q initialization for Eta for optimization

- pars:

  character vector of posterior parameters to return

- newdata:

  Default is `NULL`. If non-null, newdata is used in the uncollapse
  sampler in place of X.

- ...:

  other arguments passed to
  [pibble](https://jsilve24.github.io/fido/reference/pibble_fit.md)
  (which is used internally to fit the basset model)

- m:

  object of class bassetfit

## Value

an object of class bassetfit

## Details

the full model is given by: \$\$Y_j \sim Multinomial(\pi_j)\$\$
\$\$\pi_j = \Phi^{-1}(\eta_j)\$\$ \$\$\eta \sim MN\_{D-1 \times
N}(\Lambda, \Sigma, I_N)\$\$ \$\$\Lambda \sim GP\_{D-1 \times
Q}(\Theta(X), \Sigma, \Gamma(X))\$\$ \$\$\Sigma \sim InvWish(\upsilon,
\Xi)\$\$ Where \\\Gamma(X)\\ is short hand for the Gram matrix of the
Kernel function.

Alternatively can be used to fit an additive GP of the form: \$\$Y_j
\sim Multinomial(\pi_j)\$\$ \$\$\pi_j = \Phi^{-1}(\eta_j)\$\$ \$\$\eta
\sim MN\_{D-1 \times N}(\Lambda, \Sigma, I_N)\$\$ \$\$\Lambda =
\Lambda_1 + ... + \Lambda_p + B X\$\$ \$\$\Lambda_1 \sim GP\_{D-1 \times
Q}(\Theta_1(X), \Sigma, \Gamma_1(X))\$\$ \$\$...\$\$ \$\$\Lambda_p \sim
GP\_{D-1 \times Q}(\Theta_p(X), \Sigma, \Gamma_p(X))\$\$ \$\$B \sim
MN(\Theta_B, \Sigma, \Gamma_B)\$\$ \$\$\Sigma \sim InvWish(\upsilon,
\Xi)\$\$ Where \\\Gamma(X)\\ is short hand for the Gram matrix of the
Kernel function.

Default behavior is to use MAP estimate for uncollaping the LTP model if
laplace approximation is not preformed.

# Interface to fit orthus models

This function is largely a more user friendly wrapper around
[`optimPibbleCollapsed`](https://jsilve24.github.io/fido/reference/optimPibbleCollapsed.md)
and
[`uncollapsePibble`](https://jsilve24.github.io/fido/reference/uncollapsePibble.md)
for fitting orthus models. See details for model specification.
Notation: `N` is number of samples, `P` is the number of dimensions of
observations in the second dataset, `D` is number of multinomial
categories, `Q` is number of covariates, `iter` is the number of samples
of `eta` (e.g., the parameter `n_samples` in the function
[`optimPibbleCollapsed`](https://jsilve24.github.io/fido/reference/optimPibbleCollapsed.md))

## Usage

``` r
orthus(
  Y = NULL,
  Z = NULL,
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

- Z:

  P x N matrix of counts (if NULL uses priors only - must be
  present/absent if Y is present/absent)

- X:

  Q x N matrix of covariates (design matrix) (if NULL uses priors only,
  must be present to sample Eta)

- upsilon:

  dof for inverse wishart prior (numeric must be \> D) (default: D+3)

- Theta:

  (D-1+P) x Q matrix of prior mean for regression parameters (default:
  matrix(0, D-1+P, Q))

- Gamma:

  QxQ prior covariance matrix (default: diag(Q))

- Xi:

  (D-1+P)x(D-1+P) prior covariance matrix (default: ALR transform of
  diag(1)\*(upsilon-D)/2 - this is essentially iid on "base scale" using
  Aitchison terminology)

- init:

  (D-1) x Q initialization for Eta for optimization

- pars:

  character vector of posterior parameters to return

- ...:

  arguments passed to
  [`optimPibbleCollapsed`](https://jsilve24.github.io/fido/reference/optimPibbleCollapsed.md)
  and
  [`uncollapsePibble`](https://jsilve24.github.io/fido/reference/uncollapsePibble.md)

## Value

an object of class pibblefit

## Details

the full model is given by: \$\$Y_j \sim Multinomial(\pi_j)\$\$
\$\$\pi_j = \Phi^{-1}(\eta_j)\$\$ \$\$cbind(\eta, Z) \sim MN\_{D-1+P
\times N}(\Lambda X, \Sigma, I_N)\$\$ \$\$\Lambda \sim MN\_{D-1+P \times
Q}(\Theta, \Sigma, \Gamma)\$\$ \$\$\Sigma \sim InvWish(\upsilon,
\Xi)\$\$ Where \\\Gamma\\ is a Q x Q covariance matrix, and
\\\Phi^{-1}\\ is ALRInv_D transform. That is, the orthus model models
the latent multinomial log-ratios (Eta) and the observations of the
second dataset jointly as a linear model. This allows Sigma to also
describe the covariation between the two datasets.

Default behavior is to use MAP estimate for uncollaping the LTP model if
laplace approximation is not preformed.

## References

JD Silverman K Roche, ZC Holmes, LA David, S Mukherjee. Bayesian
Multinomial Logistic Normal Models through Marginally Latent Matrix-T
Processes. 2019, arXiv e-prints, arXiv:1903.11695

## See also

[`fido_transforms`](https://jsilve24.github.io/fido/reference/fido_transforms.md)
provide convenience methods for transforming the representation of
pibblefit objects (e.g., conversion to proportions, alr, clr, or ilr
coordinates.)

[`access_dims`](https://jsilve24.github.io/fido/reference/access_dims.md)
provides convenience methods for accessing dimensions of pibblefit
object

## Examples

``` r
sim <- orthus_sim()
fit <- orthus(sim$Y, sim$Z, sim$X)
```

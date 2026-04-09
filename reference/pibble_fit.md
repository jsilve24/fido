# Interface to fit pibble models

This function is largely a more user friendly wrapper around
[`optimPibbleCollapsed`](https://jsilve24.github.io/fido/reference/optimPibbleCollapsed.md)
and
[`uncollapsePibble`](https://jsilve24.github.io/fido/reference/uncollapsePibble.md).
See details for model specification. Notation: `N` is number of samples,
`D` is number of multinomial categories, `Q` is number of covariates,
`iter` is the number of samples of `eta` (e.g., the parameter
`n_samples` in the function
[`optimPibbleCollapsed`](https://jsilve24.github.io/fido/reference/optimPibbleCollapsed.md))

## Usage

``` r
pibble(
  Y = NULL,
  X = NULL,
  upsilon = NULL,
  Theta = NULL,
  Gamma = NULL,
  Xi = NULL,
  init = NULL,
  pars = c("Eta", "Lambda", "Sigma"),
  newdata = NULL,
  ...
)

# S3 method for class 'pibblefit'
refit(m, pars = c("Eta", "Lambda", "Sigma"), ...)
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

- newdata:

  Default is `NULL`. If non-null, newdata is used in the uncollapse
  sampler in place of X.

- ...:

  arguments passed to
  [`optimPibbleCollapsed`](https://jsilve24.github.io/fido/reference/optimPibbleCollapsed.md)
  and
  [`uncollapsePibble`](https://jsilve24.github.io/fido/reference/uncollapsePibble.md)

- m:

  object of class pibblefit

## Value

an object of class pibblefit

## Details

the full model is given by: \$\$Y_j \sim Multinomial(\pi_j)\$\$
\$\$\pi_j = \Phi^{-1}(\eta_j)\$\$ \$\$\eta \sim MN\_{D-1 \times
N}(\Lambda X, \Sigma, I_N)\$\$ \$\$\Lambda \sim MN\_{D-1 \times
Q}(\Theta, \Sigma, \Gamma)\$\$ \$\$\Sigma \sim InvWish(\upsilon,
\Xi)\$\$ Where \\\Gamma\\ is a Q x Q covariance matrix, and
\\\Phi^{-1}\\ is ALRInv_D transform.

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

Generic functions including
[`summary`](https://jsilve24.github.io/fido/reference/summary.pibblefit.md),
[`print`](https://jsilve24.github.io/fido/reference/print.pibblefit.md),
[`coef`](https://jsilve24.github.io/fido/reference/coef.pibblefit.md),
[`as.list`](https://jsilve24.github.io/fido/reference/as.list.pibblefit.md),
[`predict`](https://jsilve24.github.io/fido/reference/predict.pibblefit.md),
[`name`](https://jsilve24.github.io/fido/reference/name.pibblefit.md),
and
[`sample_prior`](https://jsilve24.github.io/fido/reference/sample_prior.pibblefit.md)
[`name_dims`](https://jsilve24.github.io/fido/reference/name_dims.md)

Plotting functions provided by
[`plot`](https://jsilve24.github.io/fido/reference/plot.pibblefit.md)
and [`ppc`](https://jsilve24.github.io/fido/reference/ppc.pibblefit.md)
(posterior predictive checks)

## Examples

``` r
sim <- pibble_sim()
fit <- pibble(sim$Y, sim$X)
```

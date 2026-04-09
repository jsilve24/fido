# Uncollapse output from optimPibbleCollapsed to full pibble Model when Sigma is known

See details for model. Should likely be called following
[`optimPibbleCollapsed`](https://jsilve24.github.io/fido/reference/optimPibbleCollapsed.md).
Notation: `N` is number of samples, `D` is number of multinomial
categories, `Q` is number of covariates, `iter` is the number of samples
of `eta` (e.g., the parameter `n_samples` in the function
`optimPibbleCollapsed`)

## Usage

``` r
uncollapsePibble_sigmaKnown(
  eta,
  X,
  Theta,
  Gamma,
  GammaComb,
  Xi,
  sigma,
  upsilon,
  seed,
  ret_mean = FALSE,
  linear = FALSE,
  ncores = -1L
)
```

## Arguments

- eta:

  array of dimension (D-1) x N x iter (e.g., `Pars` output of function
  optimPibbleCollapsed)

- X:

  matrix of covariates of dimension Q x N

- Theta:

  matrix of prior mean of dimension (D-1) x Q

- Gamma:

  covariance matrix of dimension Q x Q

- GammaComb:

  summed covariance matrix across additive components of dimension Q x
  Q.

- Xi:

  covariance matrix of dimension (D-1) x (D-1)

- sigma:

  known covariance matrix of dimension (D-1) x (D-1) x iter

- upsilon:

  scalar (must be \> D) degrees of freedom for InvWishart prior

- seed:

  seed to use for random number generation

- ret_mean:

  if true then uses posterior mean of Lambda and Sigma corresponding to
  each sample of eta rather than sampling from posterior of Lambda and
  Sigma (useful if Laplace approximation is not used (or fails) in
  optimPibbleCollapsed)

- linear:

  Boolean. Is this for a linear parameter?

- ncores:

  (default:-1) number of cores to use, if ncores==-1 then uses default
  from OpenMP typically to use all available cores.

## Value

List with components

1.  Lambda Array of dimension (D-1) x Q x iter (posterior samples)

2.  Sigma Array of dimension (D-1) x (D-1) x iter (posterior samples)

3.  The number of cores used

4.  Timer

## Details

Notation: Let \\Z_j\\ denote the J-th row of a matrix Z. While the
collapsed model is given by: \$\$Y_j \sim Multinomial(\pi_j)\$\$
\$\$\pi_j = \Phi^{-1}(\eta_j)\$\$ \$\$\eta \sim T\_{D-1, N}(\upsilon,
\Theta X, K, A)\$\$ Where \\A = I_N + X \Gamma X'\\, \\K = \Xi\\ is a
(D-1)x(D-1) covariance matrix, \\\Gamma\\ is a Q x Q covariance matrix,
and \\\Phi^{-1}\\ is ALRInv_D transform.

The uncollapsed model (Full pibble model) is given by: \$\$Y_j \sim
Multinomial(\pi_j)\$\$ \$\$\pi_j = \Phi^{-1}(\eta_j)\$\$ \$\$\eta \sim
MN\_{D-1 \times N}(\Lambda X, \Sigma, I_N)\$\$ \$\$\Lambda \sim MN\_{D-1
\times Q}(\Theta, \Sigma, \Gamma)\$\$ \$\$\Sigma \sim InvWish(\upsilon,
\Xi)\$\$ This function provides a means of sampling from the posterior
distribution of `Lambda` and `Sigma` given posterior samples of `Eta`
from the collapsed model.

## References

JD Silverman K Roche, ZC Holmes, LA David, S Mukherjee. Bayesian
Multinomial Logistic Normal Models through Marginally Latent Matrix-T
Processes. 2019, arXiv e-prints, arXiv:1903.11695

## See also

[`optimPibbleCollapsed`](https://jsilve24.github.io/fido/reference/optimPibbleCollapsed.md)

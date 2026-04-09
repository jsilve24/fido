# Function to Optimize the Collapsed Pibble Model

See details for model. Should likely be followed by function
[`uncollapsePibble`](https://jsilve24.github.io/fido/reference/uncollapsePibble.md).
Notation: `N` is number of samples, `D` is number of multinomial
categories, and `Q` is number of covariates.

## Usage

``` r
optimPibbleCollapsed(
  Y,
  upsilon,
  ThetaX,
  KInv,
  AInv,
  init,
  n_samples = 2000L,
  calcGradHess = TRUE,
  b1 = 0.9,
  b2 = 0.99,
  step_size = 0.003,
  epsilon = 1e-06,
  eps_f = 1e-10,
  eps_g = 1e-04,
  max_iter = 10000L,
  verbose = FALSE,
  verbose_rate = 10L,
  decomp_method = "cholesky",
  optim_method = "lbfgs",
  eigvalthresh = 0,
  jitter = 0,
  multDirichletBoot = -1,
  useSylv = TRUE,
  ncores = -1L,
  seed = -1L
)
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

  D-1 x D-1 precision matrix (inverse of Xi)

- AInv:

  N x N precision matrix given by \\(I_N + X'\*Gamma\*X)^{-1}\\

- init:

  D-1 x N matrix of initial guess for eta used for optimization

- n_samples:

  number of samples for Laplace Approximation (=0 very fast as no
  inversion or decomposition of Hessian is required)

- calcGradHess:

  if n_samples=0 should Gradient and Hessian still be calculated using
  closed form solutions?

- b1:

  (ADAM) 1st moment decay parameter (recommend 0.9) "aka momentum"

- b2:

  (ADAM) 2nd moment decay parameter (recommend 0.99 or 0.999)

- step_size:

  (ADAM) step size for descent (recommend 0.001-0.003)

- epsilon:

  (ADAM) parameter to avoid divide by zero

- eps_f:

  (ADAM) normalized function improvement stopping criteria

- eps_g:

  (ADAM) normalized gradient magnitude stopping criteria

- max_iter:

  (ADAM) maximum number of iterations before stopping

- verbose:

  (ADAM) if true will print stats for stopping criteria and iteration
  number

- verbose_rate:

  (ADAM) rate to print verbose stats to screen

- decomp_method:

  decomposition of hessian for Laplace approximation 'eigen' (more
  stable-slightly, slower) or 'cholesky' (less stable, faster, default)

- optim_method:

  (default:"lbfgs") or "adam"

- eigvalthresh:

  threshold for negative eigenvalues in decomposition of negative
  inverse hessian (should be \<=0)

- jitter:

  (default: 0) if \>=0 then adds that factor to diagonal of Hessian
  before decomposition (to improve matrix conditioning)

- multDirichletBoot:

  if \>0 then it overrides laplace approximation and samples eta
  efficiently at MAP estimate from pseudo Multinomial-Dirichlet
  posterior.

- useSylv:

  (default: true) if N\<D-1 uses Sylvester Determinant Identity to speed
  up calculation of log-likelihood and gradients.

- ncores:

  (default:-1) number of cores to use, if ncores==-1 then uses default
  from OpenMP typically to use all available cores.

- seed:

  (random seed for Laplace approximation – integer)

## Value

List containing (all with respect to found optima)

1.  LogLik - Log Likelihood of collapsed model (up to proportionality
    constant)

2.  Gradient - (if `calcGradHess`=true)

3.  Hessian - (if `calcGradHess`=true) of the POSITIVE LOG POSTERIOR

4.  Pars - Parameter value of eta at optima

5.  Samples - (D-1) x N x n_samples array containing posterior samples
    of eta based on Laplace approximation (if n_samples\>0)

6.  Timer - Vector of Execution Times

7.  logInvNegHessDet - the log determinant of the covariacne of the
    Laplace approximation, useful for calculating marginal likelihood

8.  logMarginalLikelihood - A calculation of the log marginal likelihood
    based on the laplace approximation

## Details

Notation: Let \\Z_j\\ denote the J-th row of a matrix Z. Model: \$\$Y_j
\sim Multinomial(\pi_j)\$\$ \$\$\pi_j = \Phi^{-1}(\eta_j)\$\$ \$\$\eta
\sim T\_{D-1, N}(\upsilon, \Theta X, K, A)\$\$ Where \\A = I_N + X
\Gamma X'\\, K is a (D-1)x(D-1) covariance matrix, \\\Gamma\\ is a Q x Q
covariance matrix, and \\\Phi^{-1}\\ is ALRInv_D transform.

Gradient and Hessian calculations are fast as they are computed using
closed form solutions. That said, the Hessian matrix can be quite large
\[N\*(D-1) x N\*(D-1)\] and storage may be an issue.

Note: Warnings about large negative eigenvalues can either signal that
the optimizer did not reach an optima or (more commonly in my
experience) that the prior / degrees of freedom for the covariance
(given by parameters `upsilon` and `KInv`) were too specific and at odds
with the observed data. If you get this warning try the following.

1.  Try restarting the optimization using a different initial guess for
    eta

2.  Try decreasing (or even increasing )`step_size` (by increments of
    0.001 or 0.002) and increasing `max_iter` parameters in optimizer.
    Also can try increasing `b1` to 0.99 and decreasing `eps_f` by a few
    orders of magnitude

3.  Try relaxing prior assumptions regarding covariance matrix. (e.g.,
    may want to consider decreasing parameter `upsilon` closer to a
    minimum value of D)

4.  Try adding small amount of jitter (e.g., set `jitter=1e-5`) to
    address potential floating point errors.

## References

S. Ruder (2016) *An overview of gradient descent optimization
algorithms*. arXiv 1609.04747

JD Silverman K Roche, ZC Holmes, LA David, S Mukherjee. *Bayesian
Multinomial Logistic Normal Models through Marginally Latent Matrix-T
Processes*. 2022, Journal of Machine Learning

## See also

[`uncollapsePibble`](https://jsilve24.github.io/fido/reference/uncollapsePibble.md)

## Examples

``` r
sim <- pibble_sim()

# Fit model for eta
fit <- optimPibbleCollapsed(sim$Y, sim$upsilon, sim$Theta%*%sim$X, sim$KInv, 
                             sim$AInv, random_pibble_init(sim$Y))  
```

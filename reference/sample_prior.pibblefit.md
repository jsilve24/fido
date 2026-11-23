# Sample from the prior distribution of pibblefit object

Note this can be used to sample from prior and then predict can be
called to get counts or LambdaX
([`predict.pibblefit`](https://jsilve24.github.io/fido/reference/predict.pibblefit.md))

## Usage

``` r
# S3 method for class 'pibblefit'
sample_prior(
  m,
  n_samples = 2000L,
  pars = c("Eta", "Lambda", "Sigma"),
  use_names = TRUE,
  ...
)
```

## Arguments

- m:

  object of class pibblefit

- n_samples:

  number of samples to produce

- pars:

  parameters to sample

- use_names:

  should names be used if available

- ...:

  currently ignored

## Value

A pibblefit object

## Details

Could be greatly speed up in the future if needed by sampling directly
from cholesky form of inverse wishart (currently implemented as header
in this library - see MatDist.h).

## Examples

``` r
# Sample prior of already fitted  pibblefit object
sim <- pibble_sim()
#> Error in chol.default(Sigma): the leading minor of order 3 is not positive
attach(sim)
#> Error: object 'sim' not found
fit <- pibble(Y, X)
#> Error: object 'X' not found
head(sample_prior(fit))
#> Error: object 'fit' not found

# Sample prior as part of model fitting
m <- pibblefit(N=as.integer(sim$N), D=as.integer(sim$D), Q=as.integer(sim$Q), 
                iter=2000L, upsilon=upsilon, 
                Xi=Xi, Gamma=Gamma, Theta=Theta, X=X, 
                coord_system="alr", alr_base=D)
#> Error: object 'sim' not found
m <- sample_prior(m)
#> Error: object 'm' not found
plot(m) # plot prior distribution (defaults to parameter Lambda) 
#> Error: object 'm' not found
```

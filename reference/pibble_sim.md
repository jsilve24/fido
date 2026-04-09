# Simulate simple pibble dataset and priors (for testing)

Simulate simple pibble dataset and priors (for testing)

## Usage

``` r
pibble_sim(D = 10, N = 30, Q = 2, use_names = TRUE, true_priors = FALSE)
```

## Arguments

- D:

  number of multinomial categories

- N:

  number of samples

- Q:

  number of covariates (first one is an intercept, must be \> 1)

- use_names:

  should samples, covariates, and categories be named

- true_priors:

  should Xi and upsilon be chosen to have mean at true simulated value

## Value

list

## Examples

``` r
sim <- pibble_sim()
```

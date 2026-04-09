# Convert pibble samples of Eta Lambda and Sigma to tidy format

Combines them all into a single tibble, see example for formatting and
column headers. Primarily designed to be used by
[`summary.pibblefit`](https://jsilve24.github.io/fido/reference/summary.pibblefit.md).

## Usage

``` r
pibble_tidy_samples(m, use_names = FALSE, as_factor = FALSE)
```

## Arguments

- m:

  an object of class pibblefit

- use_names:

  should dimension indices be replaced by dimension names if provided in
  data used to fit pibble model.

- as_factor:

  if use_names should names be returned as factor?

## Value

tibble

## Examples

``` r
sim <- pibble_sim()
fit <- pibble(sim$Y, sim$X)
fit_tidy <- pibble_tidy_samples(fit, use_names=TRUE)
head(fit_tidy)
#>   Parameter       coord sample iter        val covariate coord2
#> 1       Eta log(c1/c10)     s1    1  4.0017697      <NA>   <NA>
#> 2       Eta log(c2/c10)     s1    1  4.4819772      <NA>   <NA>
#> 3       Eta log(c3/c10)     s1    1 -4.5087071      <NA>   <NA>
#> 4       Eta log(c4/c10)     s1    1  5.5782748      <NA>   <NA>
#> 5       Eta log(c5/c10)     s1    1 -0.2482112      <NA>   <NA>
#> 6       Eta log(c6/c10)     s1    1  1.1804562      <NA>   <NA>
```

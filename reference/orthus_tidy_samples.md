# Convert orthus samples of Eta Lambda and Sigma to tidy format

Combines them all into a single tibble, see example for formatting and
column headers. Primarily designed to be used by
[`summary.orthusfit`](https://jsilve24.github.io/fido/reference/summary.orthusfit.md).

## Usage

``` r
orthus_tidy_samples(m, use_names = FALSE, as_factor = FALSE)
```

## Arguments

- m:

  an object of class orthusfit

- use_names:

  should dimension indices be replaced by dimension names if provided in
  data used to fit pibble model.

- as_factor:

  if use_names should names be returned as factor?

## Value

tibble

## Examples

``` r
sim <- orthus_sim()
fit <- orthus(sim$Y, sim$Z, sim$X)
fit_tidy <- orthus_tidy_samples(fit, use_names=TRUE)
head(fit_tidy)
#>   Parameter       coord sample iter        val covariate coord2
#> 1       Eta log(c1/c10)     s1    1  1.7949412      <NA>   <NA>
#> 2       Eta log(c2/c10)     s1    1  1.7951659      <NA>   <NA>
#> 3       Eta log(c3/c10)     s1    1 -0.5377134      <NA>   <NA>
#> 4       Eta log(c4/c10)     s1    1 -1.0676764      <NA>   <NA>
#> 5       Eta log(c5/c10)     s1    1  0.2619365      <NA>   <NA>
#> 6       Eta log(c6/c10)     s1    1 -2.0371273      <NA>   <NA>
```

# Shortcut for summarize variable with quantiles and mean

Shortcut for summarize variable with quantiles and mean

## Usage

``` r
summarise_posterior(data, var, ...)
```

## Arguments

- data:

  tidy data frame

- var:

  variable name (unquoted) to be summarised

- ...:

  other expressions to pass to summarise

## Value

data.frame

## Details

Notation: `pX` refers to the `X`% quantile

## Examples

``` r
d <- data.frame("a"=sample(1:10, 50, TRUE),
                "b"=rnorm(50))

# Summarize posterior for b over grouping of a and also calcuate
# minmum of b (in addition to normal statistics returned)
d <- dplyr::group_by(d, a)
summarise_posterior(d, b, mean.b = mean(b), min=min(b))
#> # A tibble: 10 × 9
#>        a    p2.5     p25     p50     mean     p75  p97.5   mean.b     min
#>    <int>   <dbl>   <dbl>   <dbl>    <dbl>   <dbl>  <dbl>    <dbl>   <dbl>
#>  1     1 -1.50   -1.38   -0.644  -0.615    0.123  0.316  -0.615   -1.51  
#>  2     2 -0.609  -0.311   0.0189  0.125    0.508  0.949   0.125   -0.642 
#>  3     3 -0.227  -0.0597  0.474   0.364    0.765  0.840   0.364   -0.238 
#>  4     4 -0.995  -0.306   0.714   0.290    0.904  1.12    0.290   -1.01  
#>  5     5 -1.51   -0.202   0.512  -0.00225  0.665  0.668  -0.00225 -1.65  
#>  6     6  0.483   0.496   0.591   0.751    0.846  1.29    0.751    0.482 
#>  7     7 -0.410  -0.361  -0.230  -0.0870   0.247  0.316  -0.0870  -0.415 
#>  8     8 -1.82   -0.494   0.340   0.0668   0.681  1.60    0.0668  -1.99  
#>  9     9 -0.0325  0.0389  0.284   0.761    1.13   2.27    0.761   -0.0404
#> 10    10 -0.928  -0.350  -0.106  -0.285   -0.0410 0.0535 -0.285   -0.992 
```

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
#>        a   p2.5     p25      p50    mean     p75 p97.5  mean.b    min
#>    <int>  <dbl>   <dbl>    <dbl>   <dbl>   <dbl> <dbl>   <dbl>  <dbl>
#>  1     1 -1.75  -1.24   -0.739   -0.518  -0.0201 1.09  -0.518  -1.81 
#>  2     2 -1.07  -0.0858  0.147    0.0402  0.411  0.859  0.0402 -1.18 
#>  3     3 -2.19  -1.23   -0.133   -0.0941  1.51   1.66  -0.0941 -2.29 
#>  4     4 -1.99  -1.28   -1.08    -1.000  -0.870  0.293 -1.000  -2.09 
#>  5     5 -1.72  -0.406   0.00389 -0.110   0.505  0.915 -0.110  -1.91 
#>  6     6 -1.40  -0.219  -0.165   -0.187   0.302  0.539 -0.187  -1.56 
#>  7     7  0.875  0.875   0.875    0.875   0.875  0.875  0.875   0.875
#>  8     8 -2.00  -1.12   -0.142   -0.482   0.324  0.743 -0.482  -2.09 
#>  9     9 -1.74  -0.385   0.213    0.405   1.27   2.64   0.405  -1.92 
#> 10    10  0.107  0.119   0.368    0.742   0.761  2.20   0.742   0.106
```

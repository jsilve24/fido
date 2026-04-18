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
#>        a    p2.5    p25      p50   mean    p75 p97.5 mean.b     min
#>    <int>   <dbl>  <dbl>    <dbl>  <dbl>  <dbl> <dbl>  <dbl>   <dbl>
#>  1     1 -1.74   -1.39  -0.997   -0.721 -0.190 0.537 -0.721 -1.78  
#>  2     2 -0.529   0.351  1.33     0.850  1.59  1.82   0.850 -0.627 
#>  3     3  0.0779  0.198  0.442    0.581  0.825 1.32   0.581  0.0645
#>  4     4 -0.570  -0.232  0.145    0.145  0.521 0.860  0.145 -0.608 
#>  5     5 -0.736   0.213  0.483    0.526  0.804 1.85   0.526 -0.845 
#>  6     6 -1.72   -1.14  -0.578   -0.564 -0.301 0.864 -0.564 -1.79  
#>  7     7 -1.02   -0.778 -0.00534 -0.111  0.528 0.696 -0.111 -1.03  
#>  8     8  0.243   0.680  1.17     1.27   1.81  2.38   1.27   0.195 
#>  9     9 -1.56   -0.514  0.0854  -0.114  0.302 1.14  -0.114 -1.68  
#> 10    10 -0.663  -0.424  0.433    0.383  0.719 1.94   0.383 -0.691 
```

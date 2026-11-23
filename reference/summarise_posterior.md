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
#>        a   p2.5    p25       p50    mean     p75  p97.5  mean.b    min
#>    <int>  <dbl>  <dbl>     <dbl>   <dbl>   <dbl>  <dbl>   <dbl>  <dbl>
#>  1     1 -1.17  -0.729 -0.241    -0.251  -0.103   0.921 -0.251  -1.22 
#>  2     2 -0.753 -0.703 -0.575    -0.338  -0.0686  0.366 -0.338  -0.759
#>  3     3 -0.941 -0.392 -0.111    -0.123  -0.0870  0.872 -0.123  -1.00 
#>  4     4 -0.947 -0.757 -0.000817 -0.117   0.462   0.626 -0.117  -0.952
#>  5     5 -0.194  0.304  0.563     0.557   0.816   1.30   0.557  -0.249
#>  6     6 -0.228 -0.228 -0.228    -0.228  -0.228  -0.228 -0.228  -0.228
#>  7     7 -0.698 -0.462  0.226     0.427   1.45    1.67   0.427  -0.726
#>  8     8 -0.155  0.192  0.233     0.334   0.451   0.897  0.334  -0.237
#>  9     9 -0.438 -0.399 -0.355     0.0408  0.283   0.856  0.0408 -0.443
#> 10    10 -2.66  -0.384 -0.112    -0.349   0.563   0.788 -0.349  -2.98 
```

# Print dimensions and coordinate system information for pibblefit object.

Print dimensions and coordinate system information for pibblefit object.

## Usage

``` r
# S3 method for class 'pibblefit'
print(x, summary = FALSE, ...)
```

## Arguments

- x:

  an object of class pibblefit

- summary:

  if true also calculates and prints summary

- ...:

  other arguments to pass to summary function

## Value

No direct return, prints out summary

## See also

[`summary.pibblefit`](https://jsilve24.github.io/fido/reference/summary.pibblefit.md)
summarizes posterior intervals

## Examples

``` r
sim <- pibble_sim()
fit <- pibble(sim$Y, sim$X)
print(fit)
#> pibblefit Object: 
#>   Number of Samples:      30 
#>   Number of Categories:       10 
#>   Number of Covariates:       2 
#>   Number of Posterior Samples:    2000 
#>   Contains Samples of Parameters:Eta  Lambda  Sigma
#>   Coordinate System:      alr, reference category: 10 [c10] 
#>   Log Marginal Likelihood:    -1830.873 
```

# Print dimensions and coordinate system information for orthusfit object.

Print dimensions and coordinate system information for orthusfit object.

## Usage

``` r
# S3 method for class 'orthusfit'
print(x, summary = FALSE, ...)
```

## Arguments

- x:

  an object of class orthusfit

- summary:

  if true also calculates and prints summary

- ...:

  other arguments to pass to summary function

## Value

No direct return, prints out summary

## See also

[`summary.orthusfit`](https://jsilve24.github.io/fido/reference/summary.orthusfit.md)
summarizes posterior intervals

## Examples

``` r
sim <- orthus_sim()
fit <- orthus(sim$Y, sim$Z, sim$X)
print(fit)
#> orthusfit Object: 
#>   Number of Samples:      30 
#>   Number of Categories:       10 
#>   Number of Zdimensions:  10 
#>   Number of Covariates:       2 
#>   Number of Posterior Samples:    2000 
#>   Contains Samples of Parameters:Eta  Lambda  Sigma
#>   Coordinate System:      alr, reference category: 10 [c10] 
```

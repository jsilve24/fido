# Transform Lambda into IQLR (Inter-Quantile Log-Ratio)

Takes idea from Wu et al. (citation below) and calculates IQLR for
Lambda, potentially useful if you believe there is an invariant group of
categories (e.g., taxa / genes) that are not changing (in absolute
abundance) between samples. IQLR is defined as \$\$IQLR_x =
log(x_i/g(IQVF))\$\$ for i in 1,...,D. IQVF are the CLR coordinates
whose variance is within the inter-quantile range (defined by `probs`
argument to this function). A different IQVF is fit for each posteior
sample as the IQVFs are calculted based on posterior estimates for
Lambda. The variance of a CLR coordinate is defined as the norm of each
row of Lambda\[,focus.cov\] (i.e., the covariation in Eta, explained by
those covariates). This definition of variance allows uses to exclude
variation from technical / trivial sources in calculation of IQVF/IQLR.

## Usage

``` r
lambda_to_iqlr(m, focus.cov = NULL, probs = c(0.25, 0.75))
```

## Arguments

- m:

  object of class pibblefit (e.g., output of
  [`pibble`](https://jsilve24.github.io/fido/reference/pibble_fit.md))

- focus.cov:

  vector of integers or characters specifying columns (covariates) of
  Lambda to include in calculating IQLR (if NULL, default, then uses all
  covariates)

- probs:

  bounds for categories (i.e., features / genes / taxa) to include in
  calculation of iqlr (smaller bounds means more stringent inclusion
  criteria)

## Value

array of dimension (D, Q, iter) where D is number of taxa, Q is number
of covariates, and iter is number of posterior samples.

## Details

Primarily intended for doing differential expression analysis under
assumption that only small group of categories (e.g., taxa / genes) are
changing

## References

Jia R. Wu, Jean M. Macklaim, Briana L. Genge, Gregory B. Gloor (2017)
Finding the center: corrections for asymmetry in high-throughput
sequencing datasets. arxiv:1704.01841v1

## Examples

``` r

sim <- pibble_sim()
fit <- pibble(sim$Y, sim$X)
# Use first two covariates to define iqlr, just show first 5 samples
lambda_to_iqlr(fit, 1:2)[,,1:5] 
#> , , 1
#> 
#>              [,1]        [,2]
#>  [1,]  2.66889352  4.96094014
#>  [2,] -2.71868309 -2.44737052
#>  [3,]  6.33013698  0.43797817
#>  [4,]  1.40064152 -1.02343028
#>  [5,]  0.09265135  0.43077643
#>  [6,]  0.74526606 -1.24610124
#>  [7,] -0.23411787  0.04982587
#>  [8,]  0.73380958  1.90743549
#>  [9,] -2.87971715  0.36209603
#> [10,]  1.01447965 -0.14305987
#> 
#> , , 2
#> 
#>             [,1]       [,2]
#>  [1,]  3.5110221  5.7977386
#>  [2,] -2.8121144 -2.3798702
#>  [3,]  5.5675552 -0.4590939
#>  [4,]  1.5269688 -0.6953473
#>  [5,]  0.4687443  0.0969331
#>  [6,]  1.1341357 -1.3677793
#>  [7,] -0.1014774  0.0943728
#>  [8,] -0.1403193  2.4172359
#>  [9,] -2.5207852 -0.3541093
#> [10,]  1.3174585 -0.1991465
#> 
#> , , 3
#> 
#>             [,1]        [,2]
#>  [1,]  2.8555421  5.74377741
#>  [2,] -2.7099009 -2.64106887
#>  [3,]  6.4820833  0.16940804
#>  [4,]  0.3717666 -1.25768616
#>  [5,]  0.3567992  0.76245391
#>  [6,]  1.2431116 -1.32332141
#>  [7,]  0.1268758  0.25878682
#>  [8,]  1.0685574  2.39454209
#>  [9,] -2.6834355  0.18646548
#> [10,]  1.2986605 -0.09379751
#> 
#> , , 4
#> 
#>              [,1]          [,2]
#>  [1,]  3.38108987  5.1066864494
#>  [2,] -3.48220446 -2.6604142628
#>  [3,]  6.19373102 -0.6202731907
#>  [4,]  1.49363026 -0.7799788511
#>  [5,]  0.23364431  0.3648224778
#>  [6,]  1.14710399 -0.9231080584
#>  [7,]  0.03874099  0.1659560023
#>  [8,]  0.21807167  1.9338281770
#>  [9,] -2.85880591 -0.2307412674
#> [10,]  1.49168321  0.0006412557
#> 
#> , , 5
#> 
#>              [,1]       [,2]
#>  [1,]  2.42694427  5.2836820
#>  [2,] -3.10325382 -3.2266912
#>  [3,]  6.34534913 -0.4611319
#>  [4,]  1.10872262 -1.4431521
#>  [5,]  0.62611017  0.3076211
#>  [6,]  1.24138468 -1.4403850
#>  [7,] -0.03812482 -0.1238459
#>  [8,] -0.86379483  2.2582900
#>  [9,] -1.48631247  0.6252471
#> [10,]  1.37786176 -0.5654763
#> 

```

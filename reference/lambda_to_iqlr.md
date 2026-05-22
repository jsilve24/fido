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
#>  [1,]  2.58957236  4.92307856
#>  [2,] -2.72867180 -2.43396908
#>  [3,]  6.23125673  0.39951221
#>  [4,]  1.34535203 -1.00666844
#>  [5,]  0.06578187  0.38183049
#>  [6,]  0.72442514 -1.21751979
#>  [7,] -0.32147603 -0.01583317
#>  [8,]  0.70158414  1.84930217
#>  [9,] -2.77136132  0.37488606
#> [10,]  0.96400260 -0.15262003
#> 
#> , , 2
#> 
#>             [,1]        [,2]
#>  [1,]  3.3768332  5.80569472
#>  [2,] -2.7560126 -2.35417024
#>  [3,]  5.4918800 -0.47962861
#>  [4,]  1.4832382 -0.68398137
#>  [5,]  0.4377085  0.08989149
#>  [6,]  1.0675240 -1.32505067
#>  [7,] -0.1313005  0.08077952
#>  [8,] -0.1224989  2.34840452
#>  [9,] -2.4282634 -0.33937248
#> [10,]  1.2520456 -0.18865302
#> 
#> , , 3
#> 
#>              [,1]        [,2]
#>  [1,]  2.73296991  5.78113751
#>  [2,] -2.73778260 -2.60731962
#>  [3,]  6.42522465  0.23275320
#>  [4,]  0.36566170 -1.22049490
#>  [5,]  0.34760188  0.76167290
#>  [6,]  1.16758422 -1.29842955
#>  [7,]  0.07015765  0.25397638
#>  [8,]  1.04192697  2.36040024
#>  [9,] -2.57517289  0.15852421
#> [10,]  1.24180260 -0.07277298
#> 
#> , , 4
#> 
#>              [,1]        [,2]
#>  [1,]  3.30409255  5.12908195
#>  [2,] -3.52750866 -2.68513885
#>  [3,]  6.11224988 -0.63556962
#>  [4,]  1.46781250 -0.75192820
#>  [5,]  0.20056528  0.35772761
#>  [6,]  1.08850298 -0.92166462
#>  [7,] -0.02200505  0.13463889
#>  [8,]  0.19471533  1.90166405
#>  [9,] -2.75103081 -0.22807122
#> [10,]  1.41461090 -0.01772869
#> 
#> , , 5
#> 
#>             [,1]       [,2]
#>  [1,]  2.3095180  5.3095265
#>  [2,] -3.0617301 -3.1718666
#>  [3,]  6.2833448 -0.4449522
#>  [4,]  1.0954635 -1.3571806
#>  [5,]  0.5714357  0.2982064
#>  [6,]  1.1602535 -1.4180227
#>  [7,] -0.0800954 -0.1368336
#>  [8,] -0.8042038  2.1672763
#>  [9,] -1.4515131  0.6079269
#> [10,]  1.3095067 -0.5488562
#> 

```

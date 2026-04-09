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
#>              [,1]       [,2]
#>  [1,]  2.70162174  5.7750402
#>  [2,] -2.85515471 -2.8244017
#>  [3,]  6.19878028 -1.0230645
#>  [4,]  1.02127173 -1.1836089
#>  [5,]  0.01617968  0.6125065
#>  [6,]  0.55860630 -1.4854047
#>  [7,] -0.36581855 -0.2036578
#>  [8,]  0.98956549  2.1818278
#>  [9,] -2.56944352  0.4871858
#> [10,]  0.99041411 -0.3836451
#> 
#> , , 2
#> 
#>             [,1]        [,2]
#>  [1,]  3.0112420  6.13677039
#>  [2,] -2.3726961 -2.58573845
#>  [3,]  6.5220457 -0.30854515
#>  [4,]  1.2665423 -1.72363851
#>  [5,]  0.7347088  0.46208129
#>  [6,]  1.1584673 -1.21801438
#>  [7,]  0.3014077 -0.01725855
#>  [8,]  0.5312899  3.37621095
#>  [9,] -2.9562995 -0.43455807
#> [10,]  1.5065131 -0.15119485
#> 
#> , , 3
#> 
#>             [,1]       [,2]
#>  [1,]  2.5198296  5.5127530
#>  [2,] -2.9115944 -2.2461001
#>  [3,]  7.2920404 -0.4013465
#>  [4,]  1.0641658 -0.4161133
#>  [5,]  0.3098218  0.6428295
#>  [6,]  1.1089828 -1.2640672
#>  [7,]  0.1711856  0.2486960
#>  [8,]  0.5417840  2.6207731
#>  [9,] -2.7149326 -0.9405926
#> [10,]  1.2783227 -0.2655056
#> 
#> , , 4
#> 
#>              [,1]        [,2]
#>  [1,]  3.51189478  5.19538157
#>  [2,] -3.01440896 -3.20859050
#>  [3,]  6.30172955 -0.48854871
#>  [4,]  1.08874052 -0.41983986
#>  [5,]  0.80645001  0.48683949
#>  [6,]  1.35044861 -1.09536827
#>  [7,]  0.04884086 -0.04597861
#>  [8,]  0.57179243  2.28897290
#>  [9,] -1.97108190 -1.14762602
#> [10,]  1.60116110 -0.12585754
#> 
#> , , 5
#> 
#>             [,1]       [,2]
#>  [1,]  2.7119724  5.9372748
#>  [2,] -3.2872264 -3.5465034
#>  [3,]  5.0852485 -0.1107854
#>  [4,]  0.9263608 -2.2449587
#>  [5,] -0.1595438  0.5437542
#>  [6,]  0.9076671 -1.3440209
#>  [7,] -0.4560797 -0.3794584
#>  [8,]  0.8093899  2.5817427
#>  [9,] -2.6434178  1.0072369
#> [10,]  1.0112295 -0.2445637
#> 

```

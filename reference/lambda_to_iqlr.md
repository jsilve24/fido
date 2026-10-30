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
#>  [1,]  2.67169302  4.96658199
#>  [2,] -2.71656779 -2.45115653
#>  [3,]  6.34001079  0.44269927
#>  [4,]  1.40644565 -1.02477473
#>  [5,]  0.09680193  0.42799172
#>  [6,]  0.74758121 -1.24969652
#>  [7,] -0.22939746  0.04801433
#>  [8,]  0.73698249  1.92117926
#>  [9,] -2.89100935  0.35329199
#> [10,]  1.02007607 -0.14416779
#> 
#> , , 2
#> 
#>              [,1]        [,2]
#>  [1,]  3.51868113  5.80618353
#>  [2,] -2.81439915 -2.38596493
#>  [3,]  5.57264258 -0.45371568
#>  [4,]  1.53663442 -0.69122992
#>  [5,]  0.47198303  0.09817594
#>  [6,]  1.14271789 -1.36678853
#>  [7,] -0.09644223  0.09605032
#>  [8,] -0.13913345  2.43566414
#>  [9,] -2.54021886 -0.37764569
#> [10,]  1.32410463 -0.19698736
#> 
#> , , 3
#> 
#>             [,1]        [,2]
#>  [1,]  2.8619508  5.75098271
#>  [2,] -2.7110619 -2.64668085
#>  [3,]  6.4780501  0.16421965
#>  [4,]  0.3827968 -1.25558787
#>  [5,]  0.3639917  0.76491976
#>  [6,]  1.2517695 -1.32201129
#>  [7,]  0.1330900  0.25789928
#>  [8,]  1.0696267  2.40234351
#>  [9,] -2.7041930  0.17525565
#> [10,]  1.3074184 -0.09263484
#> 
#> , , 4
#> 
#>              [,1]         [,2]
#>  [1,]  3.38926972  5.115187514
#>  [2,] -3.50230457 -2.674651863
#>  [3,]  6.20546252 -0.613841866
#>  [4,]  1.50495759 -0.774836942
#>  [5,]  0.24086423  0.366310000
#>  [6,]  1.15547016 -0.921326456
#>  [7,]  0.04815033  0.170126568
#>  [8,]  0.22203046  1.950363400
#>  [9,] -2.88245821 -0.254200002
#> [10,]  1.50148655  0.004747543
#> 
#> , , 5
#> 
#>              [,1]       [,2]
#>  [1,]  2.43174553  5.2924913
#>  [2,] -3.10632872 -3.2361631
#>  [3,]  6.35227912 -0.4664180
#>  [4,]  1.11902229 -1.4437652
#>  [5,]  0.63121390  0.3051231
#>  [6,]  1.25012465 -1.4426704
#>  [7,] -0.03071845 -0.1270803
#>  [8,] -0.87434536  2.2704608
#>  [9,] -1.49480158  0.6159748
#> [10,]  1.38580396 -0.5671379
#> 

```

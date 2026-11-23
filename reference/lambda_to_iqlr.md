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
#>  [1,]  2.67091834  4.96470801
#>  [2,] -2.72003276 -2.44916781
#>  [3,]  6.33280609  0.44027275
#>  [4,]  1.40004552 -1.02491369
#>  [5,]  0.09414965  0.43314402
#>  [6,]  0.74646241 -1.24587491
#>  [7,] -0.23266247  0.05063432
#>  [8,]  0.73511784  1.91342632
#>  [9,] -2.88162577  0.35736228
#> [10,]  1.01644191 -0.14258764
#> 
#> , , 2
#> 
#>              [,1]        [,2]
#>  [1,]  3.51342106  5.80450264
#>  [2,] -2.81187378 -2.38433824
#>  [3,]  5.56926677 -0.45542179
#>  [4,]  1.52940130 -0.69537300
#>  [5,]  0.47002213  0.09884849
#>  [6,]  1.13825196 -1.36619297
#>  [7,] -0.09909728  0.09579055
#>  [8,] -0.13902476  2.42668315
#>  [9,] -2.52862850 -0.36511718
#> [10,]  1.32035072 -0.19800092
#> 
#> , , 3
#> 
#>             [,1]        [,2]
#>  [1,]  2.8564306  5.75083913
#>  [2,] -2.7112666 -2.64321925
#>  [3,]  6.4795621  0.16797076
#>  [4,]  0.3745572 -1.25822235
#>  [5,]  0.3587575  0.76539005
#>  [6,]  1.2470281 -1.32126622
#>  [7,]  0.1293510  0.25981256
#>  [8,]  1.0699830  2.39870607
#>  [9,] -2.6915684  0.18078250
#> [10,]  1.3025359 -0.09254894
#> 
#> , , 4
#> 
#>              [,1]        [,2]
#>  [1,]  3.38477391  5.11183496
#>  [2,] -3.49605219 -2.67099842
#>  [3,]  6.20027276 -0.61729024
#>  [4,]  1.49834956 -0.77863192
#>  [5,]  0.23730900  0.36493654
#>  [6,]  1.15176288 -0.92260867
#>  [7,]  0.04390586  0.16770871
#>  [8,]  0.22004206  1.94434156
#>  [9,] -2.87015450 -0.24310097
#> [10,]  1.49699163  0.00162625
#> 
#> , , 5
#> 
#>              [,1]       [,2]
#>  [1,]  2.42767241  5.2905499
#>  [2,] -3.10353946 -3.2315830
#>  [3,]  6.34806277 -0.4633447
#>  [4,]  1.11084217 -1.4450514
#>  [5,]  0.62842707  0.3073420
#>  [6,]  1.24492090 -1.4410179
#>  [7,] -0.03463587 -0.1247014
#>  [8,] -0.86607379  2.2655170
#>  [9,] -1.48968928  0.6205523
#> [10,]  1.38102996 -0.5664318
#> 

```

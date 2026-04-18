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
#>  [1,]  2.71146763  5.7759995
#>  [2,] -2.83600181 -2.8201434
#>  [3,]  6.21761264 -1.0203479
#>  [4,]  1.03980272 -1.1813194
#>  [5,]  0.03148194  0.6107504
#>  [6,]  0.56390175 -1.4876453
#>  [7,] -0.35231497 -0.2046159
#>  [8,]  0.99318550  2.1900294
#>  [9,] -2.59688997  0.4789353
#> [10,]  1.00490204 -0.3834659
#> 
#> , , 2
#> 
#>             [,1]         [,2]
#>  [1,]  3.0292163  6.150526827
#>  [2,] -2.3634092 -2.582310950
#>  [3,]  6.5441694 -0.303325513
#>  [4,]  1.2922373 -1.710284823
#>  [5,]  0.7504182  0.466272166
#>  [6,]  1.1748118 -1.211781887
#>  [7,]  0.3199885 -0.009770112
#>  [8,]  0.5421854  3.402663821
#>  [9,] -3.0092345 -0.480597111
#> [10,]  1.5254915 -0.143235239
#> 
#> , , 3
#> 
#>             [,1]       [,2]
#>  [1,]  2.5347822  5.5245149
#>  [2,] -2.9019233 -2.2496116
#>  [3,]  7.2971256 -0.3872665
#>  [4,]  1.0844075 -0.4086840
#>  [5,]  0.3311034  0.6531207
#>  [6,]  1.1255328 -1.2612509
#>  [7,]  0.1891496  0.2564228
#>  [8,]  0.5524670  2.6479555
#>  [9,] -2.7624073 -0.9780206
#> [10,]  1.2984684 -0.2569998
#> 
#> , , 4
#> 
#>              [,1]        [,2]
#>  [1,]  3.52519216  5.21034519
#>  [2,] -3.02047160 -3.22797156
#>  [3,]  6.30534034 -0.46164556
#>  [4,]  1.11320938 -0.40412392
#>  [5,]  0.81873680  0.49524828
#>  [6,]  1.37482451 -1.09045711
#>  [7,]  0.06283319 -0.03024653
#>  [8,]  0.57166711  2.32272699
#>  [9,] -2.00932480 -1.20202334
#> [10,]  1.61770359 -0.11114129
#> 
#> , , 5
#> 
#>             [,1]       [,2]
#>  [1,]  2.7267718  5.9435322
#>  [2,] -3.2860444 -3.5519548
#>  [3,]  5.0967298 -0.1122456
#>  [4,]  0.9504458 -2.2417675
#>  [5,] -0.1428236  0.5444726
#>  [6,]  0.9229196 -1.3428317
#>  [7,] -0.4404871 -0.3804046
#>  [8,]  0.8095281  2.5917576
#>  [9,] -2.6828935  0.9928416
#> [10,]  1.0290154 -0.2432622
#> 

```

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
#>             [,1]        [,2]
#>  [1,]  2.6755143  4.96878248
#>  [2,] -2.6970678 -2.44663257
#>  [3,]  6.3564156  0.44829251
#>  [4,]  1.4198751 -1.02272623
#>  [5,]  0.1034201  0.42228032
#>  [6,]  0.7511363 -1.25486588
#>  [7,] -0.2196869  0.04522016
#>  [8,]  0.7470598  1.91874594
#>  [9,] -2.9180712  0.35884617
#> [10,]  1.0279514 -0.14515999
#> 
#> , , 2
#> 
#>              [,1]       [,2]
#>  [1,]  3.52391282  5.8089702
#>  [2,] -2.80400533 -2.3739401
#>  [3,]  5.57568396 -0.4491667
#>  [4,]  1.54761996 -0.6811336
#>  [5,]  0.47702337  0.1029327
#>  [6,]  1.15383940 -1.3634498
#>  [7,] -0.08873496  0.1023391
#>  [8,] -0.13432757  2.4304269
#>  [9,] -2.56713179 -0.3858435
#> [10,]  1.33226689 -0.1912985
#> 
#> , , 3
#> 
#>             [,1]        [,2]
#>  [1,]  2.8747825  5.74790067
#>  [2,] -2.6900039 -2.63390880
#>  [3,]  6.4720066  0.15469200
#>  [4,]  0.3991163 -1.24904759
#>  [5,]  0.3746667  0.76928431
#>  [6,]  1.2645004 -1.31726407
#>  [7,]  0.1430736  0.26001458
#>  [8,]  1.0722121  2.39624133
#>  [9,] -2.7358288  0.17007033
#> [10,]  1.3175568 -0.08904947
#> 
#> , , 4
#> 
#>              [,1]        [,2]
#>  [1,]  3.39791408  5.11780707
#>  [2,] -3.47486945 -2.65905401
#>  [3,]  6.21338910 -0.60996840
#>  [4,]  1.51781081 -0.76665130
#>  [5,]  0.25021036  0.37043177
#>  [6,]  1.16565973 -0.91475096
#>  [7,]  0.05701428  0.17591408
#>  [8,]  0.22923991  1.94832995
#>  [9,] -2.91271045 -0.26692769
#> [10,]  1.51097018  0.01186599
#> 
#> , , 5
#> 
#>              [,1]       [,2]
#>  [1,]  2.44166288  5.2909596
#>  [2,] -3.09679218 -3.2253308
#>  [3,]  6.35387959 -0.4710350
#>  [4,]  1.13432104 -1.4380439
#>  [5,]  0.63709670  0.3083269
#>  [6,]  1.26080750 -1.4393899
#>  [7,] -0.02645392 -0.1267388
#>  [8,] -0.88278213  2.2591570
#>  [9,] -1.51234641  0.6182769
#> [10,]  1.39378340 -0.5641926
#> 

```

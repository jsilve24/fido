# Return regression coefficients of orthus object

Orthus: Returned as array of dimension (D-1+P) x Q x iter (if in ALR or
ILR) otherwise (D+P) x Q x iter.

## Usage

``` r
# S3 method for class 'orthusfit'
coef(object, ...)
```

## Arguments

- object:

  an object of class orthusfit

- ...:

  other options passed to coef.orthusfit (see details)

## Value

Array of dimension (D-1) x Q x iter

## Details

Other arguments:

- use_names if column and row names were passed for Y and X in call to
  [`pibble`](https://jsilve24.github.io/fido/reference/pibble_fit.md),
  should these names be applied to output array.

# Return regression coefficients of pibblefit object

Pibble: Returned as array of dimension (D-1) x Q x iter (if in ALR or
ILR) otherwise DxQxiter (if in proportions or clr).

## Usage

``` r
# S3 method for class 'pibblefit'
coef(object, ...)
```

## Arguments

- object:

  an object of class pibblefit

- ...:

  other options passed to coef.pibblefit (see details)

## Value

Array of dimension (D-1) x Q x iter

## Details

Other arguments:

- \`use_names\` if column and row names were passed for Y and X in call
  to
  [`pibble`](https://jsilve24.github.io/fido/reference/pibble_fit.md),
  should these names be applied to output array.

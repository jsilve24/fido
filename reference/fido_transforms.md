# Transform Fit fido Parameters to other representations

These are a collection of convenience functions for transforming fido
fit objects to a number of different representations including ILR
bases, CLR coordinates, ALR coordinates, and proportions.

## Usage

``` r
to_proportions(m)

to_alr(m, d)

to_ilr(m, V = NULL)

to_clr(m)

# S3 method for class 'pibblefit'
to_proportions(m)

# S3 method for class 'orthusfit'
to_proportions(m)

# S3 method for class 'pibblefit'
to_alr(m, d)

# S3 method for class 'orthusfit'
to_alr(m, d)

# S3 method for class 'pibblefit'
to_ilr(m, V = NULL)

# S3 method for class 'orthusfit'
to_ilr(m, V = NULL)

# S3 method for class 'pibblefit'
to_clr(m)

# S3 method for class 'orthusfit'
to_clr(m)
```

## Arguments

- m:

  object of class pibblefit or orthusfit (e.g., output of
  [`pibble`](https://jsilve24.github.io/fido/reference/pibble_fit.md) or
  [`orthus`](https://jsilve24.github.io/fido/reference/orthus_fit.md))

- d:

  (integer) multinomial category to take as new alr reference

- V:

  (matrix) contrast matrix for ILR basis to transform into to (defaults
  to `create_default_ilr_base(D)`)

## Value

object

## Details

For orthus, transforms only appleid to log-ratio parameters

Note: that there is a degeneracy of representations for a covariance
matrix represented in terms of proportions. As such the function
`to_proportions` does not attempt to transform parameters Sigma or prior
Xi and instead just removes them from the pibblefit object returned.

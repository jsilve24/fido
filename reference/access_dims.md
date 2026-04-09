# Generic method for accessing model fit dimensions

Generic method for accessing model fit dimensions

## Usage

``` r
# S3 method for class 'pibblefit'
ncategories(m)

# S3 method for class 'pibblefit'
nsamples(m)

# S3 method for class 'pibblefit'
ncovariates(m)

# S3 method for class 'pibblefit'
niter(m)

# S3 method for class 'orthusfit'
ncategories(m)

# S3 method for class 'orthusfit'
nsamples(m)

# S3 method for class 'orthusfit'
ncovariates(m)

# S3 method for class 'orthusfit'
niter(m)

ncategories(m)

nsamples(m)

ncovariates(m)

niter(m)
```

## Arguments

- m:

  An object of class pibblefit

## Value

integer

## Details

An alternative approach to accessing these dimensions is to access them
directly from the pibblefit object using list indexing. \* `ncategories`
is equivalent to `m$D` \* `nsamples` is equivalent to `m$N` \*
`ncovariates` is equivalent to `m$Q`

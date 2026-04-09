# Generic method for getting and setting dimension names of fit object

Generic method for getting and setting dimension names of fit object

## Usage

``` r
# S3 method for class 'pibblefit'
names_covariates(m)

# S3 method for class 'pibblefit'
names_samples(m)

# S3 method for class 'pibblefit'
names_categories(m)

# S3 method for class 'pibblefit'
names_coords(m)

# S3 method for class 'pibblefit'
names_covariates(m) <- value

# S3 method for class 'pibblefit'
names_samples(m) <- value

# S3 method for class 'pibblefit'
names_categories(m) <- value

names_covariates(m)

names_samples(m)

names_categories(m)

names_coords(m)

names_covariates(m) <- value

names_samples(m) <- value

names_categories(m) <- value
```

## Arguments

- m:

  object

- value:

  character vector (or NULL)

## Value

A vector of names

## Details

`names_coords` is different than `names_categories`. `names_categories`
provides access to the basic names of each multinomial category. In
contrast, `names_coords` provides access to the names of the coordinates
in which an object is represented. These coordinate names are based on
the category names. For example, category names may be, (OTU1, ...,
OTUD) where as coordinate names could be (log(OTU1/OTUD), etc...) if
object is in default coordinate system.

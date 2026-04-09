# Check vector/matrix/data.frame for expected dimensions or throw error

Check vector/matrix/data.frame for expected dimensions or throw error

## Usage

``` r
check_dims(x, d, par)
```

## Arguments

- x:

  object to check

- d:

  expected dimensions

- par:

  character name of x (for error message)

## Value

nothing if no error, otherwise throws error

## Examples

``` r
y <- c(1,3,4)
check_dims(y, 3, "y")
```

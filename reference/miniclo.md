# Closure operator

Closure operator

## Usage

``` r
miniclo(x)
```

## Arguments

- x:

  vector or matrix (rows are samples, parts are columns) of data in
  simplex

## Value

x with row entries divided by sum of row (converts vectors to row
matricies)

## Examples

``` r
x <- matrix(runif(30), 10, 3)
x <- miniclo(x)
```

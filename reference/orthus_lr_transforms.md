# Log-Ratio transforms for orthus objects

Log-Ratio transforms for orthus objects

## Usage

``` r
oglr(x, s, V)

oglrInv(x, s, V)

oalr(x, s, d = NULL)

oalrInv(y, s, d = NULL)

oilr(x, s, V = NULL)

oilrInv(y, s, V = NULL)

oclr(x, s)

oclrInv(x, s)
```

## Arguments

- x:

  orthus data array (e.g., first s rows are multinomial parameters or
  log-ratios)

- s:

  first s rows of x are transformed

- V:

  transformation matrix (defines transform)

- d:

  for ALR, which component (integer position) to take as reference
  (default is ncol(x)) for alrInv corresponds to column position in
  untransformed matrix.

- y:

  orthus data array (e.g., first s rows are multinomial parameters or
  log-ratios)

## Value

A matrix

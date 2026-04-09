# Convert orthus covariance matricies between representations

Convert orthus covariance matricies between representations

## Usage

``` r
oilrvar2ilrvar(Sigma, s, V1, V2)

oilrvar2clrvar(Sigma, s, V)

oclrvar2ilrvar(Sigma, s, V)

oalrvar2clrvar(Sigma, s, d1)

oclrvar2alrvar(Sigma, s, d2)

oalrvar2alrvar(Sigma, s, d1, d2)

oalrvar2ilrvar(Sigma, s, d1, V2)

oilrvar2alrvar(Sigma, s, V1, d2)
```

## Arguments

- Sigma:

  covariance matrix arrat in specified transformed space
  (dim(Sigma)\[3\]=iter)

- s:

  first s rows and colums of Sigma are transformed

- V1:

  ILR contrast matrix of basis Sigma is already in

- V2:

  ILR contrast matrix of basis Sigma is desired in

- V:

  ILR contrast matrix (i.e., transformation matrix of ILR)

- d1:

  alr reference element Sigma is already expressed with respec to

- d2:

  alr reference element Sigma is to be expressed with respect to

## Value

matrix

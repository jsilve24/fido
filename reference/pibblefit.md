# Create pibblefit object

Create pibblefit object

## Usage

``` r
pibblefit(
  D,
  N,
  Q,
  coord_system,
  iter = NULL,
  alr_base = NULL,
  ilr_base = NULL,
  Eta = NULL,
  Lambda = NULL,
  Sigma = NULL,
  Sigma_default = NULL,
  Y = NULL,
  X = NULL,
  upsilon = NULL,
  Theta = NULL,
  Xi = NULL,
  Xi_default = NULL,
  Gamma = NULL,
  init = NULL,
  names_categories = NULL,
  names_samples = NULL,
  names_covariates = NULL
)
```

## Arguments

- D:

  number of multinomial categories

- N:

  number of samples

- Q:

  number of covariates

- coord_system:

  coordinate system objects are represented in (options include "alr",
  "clr", "ilr", and "proportions")

- iter:

  number of posterior samples

- alr_base:

  integer category used as reference (required if coord_system=="alr")

- ilr_base:

  (D x D-1) contrast matrix (required if coord_system=="ilr")

- Eta:

  Array of samples of Eta

- Lambda:

  Array of samples of Lambda

- Sigma:

  Array of samples of Sigma (null if coord_system=="proportions")

- Sigma_default:

  Array of samples of Sigma in alr base D, used if
  coord_system=="proportions"

- Y:

  DxN matrix of observed counts

- X:

  QxN design matrix

- upsilon:

  scalar prior dof of inverse wishart prior

- Theta:

  prior mean of Lambda

- Xi:

  Matrix of prior covariance for inverse wishart (null if
  coord_system=="proportions")

- Xi_default:

  Matrix of prior covariance for inverse wishart in alr base D (used if
  coord_system=="proportions")

- Gamma:

  QxQ covariance matrix prior for Lambda

- init:

  matrix initial guess for Lambda used for optimization

- names_categories:

  character vector

- names_samples:

  character vector

- names_covariates:

  character vector

## Value

object of class pibblefit

## See also

[`pibble`](https://jsilve24.github.io/fido/reference/pibble_fit.md)

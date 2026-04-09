# Predict using basset

Predict using basset

## Usage

``` r
# S3 method for class 'bassetfit'
predict(
  object,
  newdata = NULL,
  response = "Lambda",
  size = NULL,
  use_names = TRUE,
  summary = FALSE,
  iter = NULL,
  from_scratch = FALSE,
  ...
)
```

## Arguments

- object:

  An object of class pibblefit

- newdata:

  An optional matrix for which to evaluate prediction.

- response:

  Options = "Lambda":Mean of regression, "Eta", "Y": counts

- size:

  the number of counts per sample if response="Y" (as vector or matrix),
  default if newdata=NULL and response="Y" is to use colsums of m\$Y.
  Otherwise uses median colsums of object\$Y as default. If passed as a
  matrix should have dimensions ncol(newdata) x iter.

- use_names:

  if TRUE apply names to output

- summary:

  if TRUE, posterior summary of predictions are returned rather than
  samples

- iter:

  number of iterations to return if NULL uses object\$iter

- from_scratch:

  should predictions of Y come from fitted Eta or from predictions of
  Eta from posterior of Lambda? (default: false)

- ...:

  other arguments passed to summarise_posterior

## Value

(if summary==FALSE) array D x N x iter; (if summary==TRUE) tibble with
calculated posterior summaries

## Details

currently only implemented for pibblefit objects in coord_system
"default" "alr", or "ilr".

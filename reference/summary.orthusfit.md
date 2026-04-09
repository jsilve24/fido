# Summarise orthusfit object and print posterior quantiles

Default calculates median, mean, 50% and 95% credible interval

## Usage

``` r
# S3 method for class 'orthusfit'
summary(
  object,
  pars = NULL,
  use_names = TRUE,
  as_factor = FALSE,
  gather_prob = FALSE,
  ...
)
```

## Arguments

- object:

  an object of class orthusfit

- pars:

  character vector (default: c("Eta", "Lambda", "Sigma"))

- use_names:

  should summary replace dimension indices with orthusfit names if names
  Y and X were named in call to
  [`orthus`](https://jsilve24.github.io/fido/reference/orthus_fit.md)

- as_factor:

  if use_names and as_factor then returns names as factors (useful for
  maintaining orderings when plotting)

- gather_prob:

  if TRUE then prints quantiles in long format rather than wide (useful
  for some plotting functions)

- ...:

  other expressions to pass to summarise (using name 'val' unquoted is
  probably what you want)

## Value

A list

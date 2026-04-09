# Holds information on coordinates system to later be reapplied

`store_coord` stores coordinate information for pibblefit object and can
be reapplied with function `reapply_coord`. Some coordinate systems are
not useful for computation and this makes it simple keep returned object
from computations in the same coordinate system as the input. (Likely
most useful inside of a package)

## Usage

``` r
store_coord(m)

reapply_coord(m, l)
```

## Arguments

- m:

  object of class pibblefit

- l:

  object returned by function `store_coord`

## Value

`store_coord` list with important information to identify c coordinate
system of pibblefit object. `reapply_coord` pibblefit object in
coordinate system previously stored.

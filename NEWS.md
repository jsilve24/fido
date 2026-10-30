# fido 1.1.5

This release fixes a long-standing error in how the prior degrees of freedom
(`upsilon`) were used. **Results from `pibble()`, `basset()` and `orthus()`
fitted with earlier versions may change; please refit.**

## What was wrong

fido fits models in two steps. It first finds `Eta` using a "collapsed" model in
which `Lambda` and `Sigma` have been integrated out. It then draws `Lambda` and
`Sigma` given `Eta` (the "uncollapse" step). Both steps should use the same
inverse-Wishart prior on `Sigma`, with `upsilon` degrees of freedom.

The collapsed step expects its degrees of freedom on a different scale (the
matrix-t scale). Converting the inverse-Wishart `upsilon` to that scale requires
subtracting `D - 2`. Since 2018, `pibble()` passed `upsilon` to the collapsed
step without this conversion. The collapsed step therefore used a prior on
`Sigma` that was more informative than the one requested. The uncollapse step
used the correct prior, so the two steps did not agree.

## Who is affected

* `pibble()` and `basset()` (including additive `basset()` models) with three or
  more categories (`D > 2`). With `D = 2` the two scales coincide and results are
  unchanged. The larger `D` is relative to `upsilon`, the larger the effect.
* `orthus()` in all cases, because it had a similar conversion error (`upsilon + P`
  was used where `upsilon - D + 2` is correct).
* `maltipoo()` (experimental) received the same fix.

## Other fixes

* The log marginal likelihood reported by the Laplace approximation
  (`logMarginalLikelihood`) was missing several constant terms and is now
  correctly normalized. Stored values from earlier versions (for any `D`) should
  not be compared with new ones.
* `pibble()`, `basset()`, `orthus()` and `maltipoo()` now stop with an error if
  `upsilon` is too small for a valid inverse-Wishart prior (`upsilon <= D-2`,
  or `upsilon <= D+P-2` for `orthus()`). Previously such values were accepted
  silently and produced meaningless fits.
* `orthus_sim(true_priors = FALSE)` now returns `upsilon = D+P` (previously `D`,
  which is not a valid prior for the `(D-1+P)`-dimensional covariance).
* `orthus_sim(true_priors = TRUE)` now uses `upsilon = D+P+10` so the
  inverse-Wishart prior mean equals the simulated `Sigma`, including with
  the default dimensions. Previously the scale was zero or negative when
  `P >= 10`.
* `orthus_sim()` now preserves matrix dimensions with two categories or a
  single Gaussian variable, sample, or covariate. Invalid dimensions are
  rejected before simulation, and its fixed covariance entry no longer makes
  the simulated covariance singular when both corresponding variances are one.
* Non-integer `upsilon` was silently rounded down when drawing `Sigma`; it is
  now used as given.
* `ret_mean = TRUE` now errors when `upsilon + N` is too small for the
  conditional posterior mean of `Sigma` to exist. Use posterior draws
  (`ret_mean = FALSE`) with these proper, heavy-tailed posteriors.
* `basset()` with more than one component (lists of `Theta` and `Gamma`) passed
  `ncores` into the `ret_mean` slot of the component sampler. When `ncores` was
  set this could return component means instead of posterior draws, and the
  `seed` and `ncores` settings were ignored. Arguments are now passed by name,
  and each component uses its own random stream, so component uncertainty is
  now correct.
* Added `fido_build_info()` to report native runtime threading configuration and
  compiler support for OpenMP, Eigen, and MKL.
* Default thread limits are capped at two when R sets
  `_R_CHECK_LIMIT_CORES_=TRUE` for package checks.
* Added opt-in support for compiling with Intel MKL and GNU OpenMP via
  `--configure-args="--with-mkl"` on supported Linux x86_64 systems, with automatic
  thread safety guards during parallel uncollapsing.
* Fixed error handling in `LaplaceApproximation` to verify Cholesky factor status
  and preserve block-diagonal partial Hessian handling when MKL is disabled.

## What you need to do

* **Refit** analyses from the original data and priors. Old fits can still be
  loaded, but there is no way to correct their draws by adjusting `upsilon`
  after the fact; predictions from old fit objects reflect the old results.
* The meaning of the `upsilon` argument has not changed: it is the standard
  inverse-Wishart degrees of freedom. For `pibble()` and `basset()`, the prior
  mean of `Sigma` is `Xi / (upsilon - D)` when `upsilon > D`. For `orthus()`,
  it is `Xi / (upsilon - D - P)` when `upsilon > D + P`. Existing scripts need
  no edits.
* If you call the low-level functions `optimPibbleCollapsed()`,
  `loglikPibbleCollapsed()`, `gradPibbleCollapsed()` or `hessPibbleCollapsed()`
  directly, note that they take matrix-t degrees of freedom. Pass
  `upsilon - D + 2` rather than `upsilon` (see their help pages).
  `uncollapsePibble()` still takes `upsilon` unchanged.

# fido 1.1.4

* Fixed the posterior mean of `Sigma` returned by `uncollapsePibble(..., ret_mean = TRUE)` to use the correct inverse-Wishart mean.
* Fixed OpenMP sampling in the uncollapse routines so posterior draws are reproducible across `ncores` for a fixed seed and no longer depend on overlapping thread RNG streams.
* Added input validation to the multinomial-Dirichlet bootstrap so Dirichlet concentration parameters must be positive and finite.

# fido 1.1.3

* Replaced explicit integer sequences with `seq_len()` and `seq_along()` for improved robustness.

# fido 1.1.1

* Fixing a minor bug in the configure script

# fido 1.1.0

* Updating the basset model to allow for additive Gaussian processes.
* Updating the random seed generation.
* Removed maltipoo as an export. There were a few small changes to names of internal maltipoo functions (e.g., "verify.maltipoo" became "verify_maltipoo").
* Updating the configure script to a compile but that occured when OpenMP is present but R is not configured to use it.
* Fixed a few small bugs (predicting with one sample; wrong constant; etc.).

# fido 1.0.4

* changing default optimizer (lbfgs rather than adam), adam is more sensitive to random seed
* adding new r2 methods for calculating Bayesian R2 for pibble and basset models
* predict.bassetfit now accepts a newdata=NULL argument to predict on training set

bug fixes

* fixing a few discrepancies with base R commands (print, summary, etc).


# fido 1.0.3

bug fixes

* Updating documentation under the newest version of R (4.2.1)

# fido 1.0.2

bug fixes

* Removing a timing test because it is irrelevant and can fail under certain set-ups.

# fido 1.0.1

bug fixes 

* Changing phyloseq from imports to depends, was causing issues installing from CRAN. 
* Fixing bug in configure.ac that was causing openmp parallelization to not work on some machines. 

# fido 1.0.0

More changes to make the final push for CRAN! Changes include:

* Updating the description file to match CRAN's standards
* Removing the dependency to driver
* Linking the one remote package (MicrobeDS) to a separate CRAN-like Github repo using the "drat" package. This is to match CRAN's policy against no use of "Remotes".
* Precomputing the longer running vignettes, cutting R CMD check times in half. The original vignettes can be found in the vignettes folder with a ".orig" extension for reproducibility.

# fido 0.1.14

* More changes to prepare for version 0.2 (and CRAN soon), including:
* Lots of cleaning of small warnings and notes
* Fixing erorrs with the base plot function
* Cleaning up some manual entries
* Fixing the errors posted on Github (to_ilr among others)

# fido 0.1.13

* tons of tiny changes to prepare for version 0.2 (and ultimately CRAN) featured changes include:
* plot and coef methods for *fit objects now abide by CRAN generic naming scheme - may cause
  some problems if prior code had positional arguments to these functions
* now being tested using travis integration (vignettes not tested there due to time constraints)
* merged Kim's fixes to maltipoo code (will almost certainly change maltipoo results; there
  had been a bug in prior versions)
* Internally stored mallard and mallard_family data now import as `mallard` and `mallard_family` 
  rather than both importing as `ps`. Solved a R CMD CHECK warning. 
* Lots of tiny updates to documentation 

# fido 0.1.12

* orthus added for jointly fitting multinomial and gaussian data (e.g., 16S microbiome and metabolomics data)

# fido 0.1.11

* Kim fixed windows installation problems (Thanks Kim!)
* Fixed bug in predict that was effecting ppc
* Fixed error when pars!=NULL in summary
* lambda_to_iqlr now had default behavior when focus.cov=NULL
* Fixed issue with Xi=NULL in pibblefit when trying to transform
* small changes to make pibble wrapper around c++ functions faster

# fido 0.1.10

* basset added for fitting non-linear regression using fido
* numerous changes to make installation easier
* Added a `NEWS.md` file to track changes to the package.

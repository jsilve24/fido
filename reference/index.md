# Package index

## Model functions

Functions for fitting the main models

- [`pibble()`](https://jsilve24.github.io/fido/reference/pibble_fit.md)
  [`refit(`*`<pibblefit>`*`)`](https://jsilve24.github.io/fido/reference/pibble_fit.md)
  : Interface to fit pibble models
- [`refit()`](https://jsilve24.github.io/fido/reference/refit.md) :
  Generic method for fitting model from passed model fit object
- [`basset()`](https://jsilve24.github.io/fido/reference/basset_fit.md)
  [`refit(`*`<bassetfit>`*`)`](https://jsilve24.github.io/fido/reference/basset_fit.md)
  : Interface to fit basset models
- [`orthus()`](https://jsilve24.github.io/fido/reference/orthus_fit.md)
  : Interface to fit orthus models

## Model checking

- [`plot(`*`<pibblefit>`*`)`](https://jsilve24.github.io/fido/reference/plot.pibblefit.md)
  : Plot Summaries of Posterior Distribution of pibblefit Parameters
- [`ppc()`](https://jsilve24.github.io/fido/reference/ppc.md) : Generic
  method for visualizing posterior predictive checks
- [`ppc(`*`<pibblefit>`*`)`](https://jsilve24.github.io/fido/reference/ppc.pibblefit.md)
  : Visualization of Posterior Predictive Check of fit model
- [`ppc_summary()`](https://jsilve24.github.io/fido/reference/ppc_summary.md)
  : Generic Method to Plot Posterior Predictive Summaries

## Useful functions for determining model fit

- [`plot(`*`<pibblefit>`*`)`](https://jsilve24.github.io/fido/reference/plot.pibblefit.md)
  : Plot Summaries of Posterior Distribution of pibblefit Parameters
- [`ppc()`](https://jsilve24.github.io/fido/reference/ppc.md) : Generic
  method for visualizing posterior predictive checks
- [`ppc(`*`<pibblefit>`*`)`](https://jsilve24.github.io/fido/reference/ppc.pibblefit.md)
  : Visualization of Posterior Predictive Check of fit model
- [`ppc_summary()`](https://jsilve24.github.io/fido/reference/ppc_summary.md)
  : Generic Method to Plot Posterior Predictive Summaries
- [`r2()`](https://jsilve24.github.io/fido/reference/r2.md) : Generic
  Method to Calculate R2 for Fitted Model
- [`sample_prior()`](https://jsilve24.github.io/fido/reference/sample_prior.md)
  : Generic method for sampling from prior distribution of object
- [`sample_prior(`*`<pibblefit>`*`)`](https://jsilve24.github.io/fido/reference/sample_prior.pibblefit.md)
  : Sample from the prior distribution of pibblefit object
- [`summarise_posterior()`](https://jsilve24.github.io/fido/reference/summarise_posterior.md)
  : Shortcut for summarize variable with quantiles and mean

## Transformation functions

- [`to_proportions()`](https://jsilve24.github.io/fido/reference/fido_transforms.md)
  [`to_alr()`](https://jsilve24.github.io/fido/reference/fido_transforms.md)
  [`to_ilr()`](https://jsilve24.github.io/fido/reference/fido_transforms.md)
  [`to_clr()`](https://jsilve24.github.io/fido/reference/fido_transforms.md)
  : Transform Fit fido Parameters to other representations

## Helper functions

- [`alr()`](https://jsilve24.github.io/fido/reference/alr.md) : Compute
  the ALR of a matrix
- [`alrInv()`](https://jsilve24.github.io/fido/reference/alrInv.md) :
  Compute the inverse ALR of a matrix
- [`alrInv_array()`](https://jsilve24.github.io/fido/reference/alrInv_array.md)
  : Compute the ALR of an array
- [`alr_array()`](https://jsilve24.github.io/fido/reference/alr_array.md)
  : Compute the ALR of an array
- [`as.list(`*`<orthusfit>`*`)`](https://jsilve24.github.io/fido/reference/as.list.orthusfit.md)
  : Convert object of class orthusfit to a list
- [`as.list(`*`<pibblefit>`*`)`](https://jsilve24.github.io/fido/reference/as.list.pibblefit.md)
  : Convert object of class pibblefit to a list
- [`as.orthusfit()`](https://jsilve24.github.io/fido/reference/as.orthusfit.md)
  : convert list to orthusfit
- [`as.pibblefit()`](https://jsilve24.github.io/fido/reference/as.pibblefit.md)
  : convert list to pibblefit
- [`check_dims()`](https://jsilve24.github.io/fido/reference/check_dims.md)
  : Check vector/matrix/data.frame for expected dimensions or throw
  error
- [`clr_array()`](https://jsilve24.github.io/fido/reference/clr_array.md)
  : Compute the CLR of an array
- [`coef(`*`<orthusfit>`*`)`](https://jsilve24.github.io/fido/reference/coef.orthusfit.md)
  : Return regression coefficients of orthus object
- [`coef(`*`<pibblefit>`*`)`](https://jsilve24.github.io/fido/reference/coef.pibblefit.md)
  : Return regression coefficients of pibblefit object
- [`conjugateLinearModel()`](https://jsilve24.github.io/fido/reference/conjugateLinearModel.md)
  : Solve Bayesian Multivariate Conjugate Linear Model
- [`oilrvar2ilrvar()`](https://jsilve24.github.io/fido/reference/convert_orthus_covariance.md)
  [`oilrvar2clrvar()`](https://jsilve24.github.io/fido/reference/convert_orthus_covariance.md)
  [`oclrvar2ilrvar()`](https://jsilve24.github.io/fido/reference/convert_orthus_covariance.md)
  [`oalrvar2clrvar()`](https://jsilve24.github.io/fido/reference/convert_orthus_covariance.md)
  [`oclrvar2alrvar()`](https://jsilve24.github.io/fido/reference/convert_orthus_covariance.md)
  [`oalrvar2alrvar()`](https://jsilve24.github.io/fido/reference/convert_orthus_covariance.md)
  [`oalrvar2ilrvar()`](https://jsilve24.github.io/fido/reference/convert_orthus_covariance.md)
  [`oilrvar2alrvar()`](https://jsilve24.github.io/fido/reference/convert_orthus_covariance.md)
  : Convert orthus covariance matricies between representations
- [`create_default_ilr_base()`](https://jsilve24.github.io/fido/reference/create_default_ilr_base.md)
  : Create a default ILR base
- [`gather_array()`](https://jsilve24.github.io/fido/reference/gather_array.md)
  : Gather Multidimensional Array to Tidy Tibble
- [`SE()`](https://jsilve24.github.io/fido/reference/kernels.md)
  [`LINEAR()`](https://jsilve24.github.io/fido/reference/kernels.md) :
  Multivariate RBF Kernel
- [`lambda_to_iqlr()`](https://jsilve24.github.io/fido/reference/lambda_to_iqlr.md)
  : Transform Lambda into IQLR (Inter-Quantile Log-Ratio)
- [`loglikPibbleCollapsed()`](https://jsilve24.github.io/fido/reference/loglikPibbleCollapsed.md)
  [`gradPibbleCollapsed()`](https://jsilve24.github.io/fido/reference/loglikPibbleCollapsed.md)
  [`hessPibbleCollapsed()`](https://jsilve24.github.io/fido/reference/loglikPibbleCollapsed.md)
  : Calculations for the Collapsed Pibble Model
- [`lmvgamma()`](https://jsilve24.github.io/fido/reference/lmvgamma.md)
  : Log of Multivarate Gamma Function - Gamma_p(a)
- [`lmvgamma_deriv()`](https://jsilve24.github.io/fido/reference/lmvgamma_deriv.md)
  : Derivative of Log of Multivariate Gamma Function - Gamma_p(a)
- [`miniclo()`](https://jsilve24.github.io/fido/reference/miniclo.md) :
  Closure operator
- [`miniclo_array()`](https://jsilve24.github.io/fido/reference/miniclo_array.md)
  : Closure Operation applied to array on margin
- [`mongrel()`](https://jsilve24.github.io/fido/reference/mongrel-deprecated.md)
  : mongrel
- [`name()`](https://jsilve24.github.io/fido/reference/name.md) :
  Generic method for applying names to an object
- [`name(`*`<orthusfit>`*`)`](https://jsilve24.github.io/fido/reference/name.orthusfit.md)
  : S3 for orthusfit apply names to orthusfit object
- [`name(`*`<pibblefit>`*`)`](https://jsilve24.github.io/fido/reference/name.pibblefit.md)
  : S3 for pibblefit apply names to pibblefit object
- [`names_covariates()`](https://jsilve24.github.io/fido/reference/name_dims.md)
  [`names_samples()`](https://jsilve24.github.io/fido/reference/name_dims.md)
  [`names_categories()`](https://jsilve24.github.io/fido/reference/name_dims.md)
  [`names_coords()`](https://jsilve24.github.io/fido/reference/name_dims.md)
  [`` `names_covariates<-`() ``](https://jsilve24.github.io/fido/reference/name_dims.md)
  [`` `names_samples<-`() ``](https://jsilve24.github.io/fido/reference/name_dims.md)
  [`` `names_categories<-`() ``](https://jsilve24.github.io/fido/reference/name_dims.md)
  : Generic method for getting and setting dimension names of fit object
- [`ncategories()`](https://jsilve24.github.io/fido/reference/access_dims.md)
  [`nsamples()`](https://jsilve24.github.io/fido/reference/access_dims.md)
  [`ncovariates()`](https://jsilve24.github.io/fido/reference/access_dims.md)
  [`niter()`](https://jsilve24.github.io/fido/reference/access_dims.md)
  : Generic method for accessing model fit dimensions
- [`optimPibbleCollapsed()`](https://jsilve24.github.io/fido/reference/optimPibbleCollapsed.md)
  : Function to Optimize the Collapsed Pibble Model
- [`orthusfit()`](https://jsilve24.github.io/fido/reference/orthusfit.md)
  : Create orthusfit object
- [`oglr()`](https://jsilve24.github.io/fido/reference/orthus_lr_transforms.md)
  [`oglrInv()`](https://jsilve24.github.io/fido/reference/orthus_lr_transforms.md)
  [`oalr()`](https://jsilve24.github.io/fido/reference/orthus_lr_transforms.md)
  [`oalrInv()`](https://jsilve24.github.io/fido/reference/orthus_lr_transforms.md)
  [`oilr()`](https://jsilve24.github.io/fido/reference/orthus_lr_transforms.md)
  [`oilrInv()`](https://jsilve24.github.io/fido/reference/orthus_lr_transforms.md)
  [`oclr()`](https://jsilve24.github.io/fido/reference/orthus_lr_transforms.md)
  [`oclrInv()`](https://jsilve24.github.io/fido/reference/orthus_lr_transforms.md)
  : Log-Ratio transforms for orthus objects
- [`orthus_tidy_samples()`](https://jsilve24.github.io/fido/reference/orthus_tidy_samples.md)
  : Convert orthus samples of Eta Lambda and Sigma to tidy format
- [`pibblefit()`](https://jsilve24.github.io/fido/reference/pibblefit.md)
  : Create pibblefit object
- [`predict(`*`<bassetfit>`*`)`](https://jsilve24.github.io/fido/reference/predict.bassetfit.md)
  : Predict using basset
- [`predict(`*`<pibblefit>`*`)`](https://jsilve24.github.io/fido/reference/predict.pibblefit.md)
  : Predict response from new data
- [`print(`*`<orthusfit>`*`)`](https://jsilve24.github.io/fido/reference/print.orthusfit.md)
  : Print dimensions and coordinate system information for orthusfit
  object.
- [`print(`*`<pibblefit>`*`)`](https://jsilve24.github.io/fido/reference/print.pibblefit.md)
  : Print dimensions and coordinate system information for pibblefit
  object.
- [`pibble_tidy_samples()`](https://jsilve24.github.io/fido/reference/pibble_tidy_samples.md)
  : Convert pibble samples of Eta Lambda and Sigma to tidy format
- [`random_pibble_init()`](https://jsilve24.github.io/fido/reference/random_pibble_init.md)
  : Provide random initialization for pibble model
- [`req()`](https://jsilve24.github.io/fido/reference/req.md) : Generic
  method for ensuring object contains required elements
- [`req(`*`<orthusfit>`*`)`](https://jsilve24.github.io/fido/reference/req.orthusfit.md)
  : require elements to be non-null in orthusfit or throw error
- [`req(`*`<pibblefit>`*`)`](https://jsilve24.github.io/fido/reference/req.pibblefit.md)
  : require elements to be non-null in pibblefit or throw error
- [`store_coord()`](https://jsilve24.github.io/fido/reference/store_coord.md)
  [`reapply_coord()`](https://jsilve24.github.io/fido/reference/store_coord.md)
  : Holds information on coordinates system to later be reapplied
- [`summary(`*`<orthusfit>`*`)`](https://jsilve24.github.io/fido/reference/summary.orthusfit.md)
  : Summarise orthusfit object and print posterior quantiles
- [`summary(`*`<pibblefit>`*`)`](https://jsilve24.github.io/fido/reference/summary.pibblefit.md)
  : Summarise pibblefit object and print posterior quantiles
- [`uncollapsePibble()`](https://jsilve24.github.io/fido/reference/uncollapsePibble.md)
  : Uncollapse output from optimPibbleCollapsed to full pibble Model
- [`uncollapsePibble_sigmaKnown()`](https://jsilve24.github.io/fido/reference/uncollapsePibble_sigmaKnown.md)
  : Uncollapse output from optimPibbleCollapsed to full pibble Model
  when Sigma is known
- [`verify()`](https://jsilve24.github.io/fido/reference/verify.md) :
  Generic method for verifying new objects
- [`verify(`*`<bassetfit>`*`)`](https://jsilve24.github.io/fido/reference/verify.bassetfit.md)
  : Simple verification of passed bassetfit object
- [`verify(`*`<orthusfit>`*`)`](https://jsilve24.github.io/fido/reference/verify.orthusfit.md)
  : Simple verification of passed orthusfit object
- [`verify(`*`<pibblefit>`*`)`](https://jsilve24.github.io/fido/reference/verify.pibblefit.md)
  : Simple verification of passed pibblefit object

## Data sets and data simulators

- [`fido`](https://jsilve24.github.io/fido/reference/fido_package.md)
  [`fido-package`](https://jsilve24.github.io/fido/reference/fido_package.md)
  [`fido_package`](https://jsilve24.github.io/fido/reference/fido_package.md)
  : fido: Fitting and Analysis of Multinomial Logistic Normal Models
- [`mallard`](https://jsilve24.github.io/fido/reference/mallard.md) :
  Data from Silverman et al. (2018) Microbiome
- [`mallard_family`](https://jsilve24.github.io/fido/reference/mallard_family.md)
  : Data from Silverman et al. (2018) Microbiome
- [`metadata`](https://jsilve24.github.io/fido/reference/metadata.md) :
  Data from Silverman et al. (2019) bioRxiv
- [`orthus_sim()`](https://jsilve24.github.io/fido/reference/orthus_sim.md)
  : Simulate simple orthus dataset and priors (for testing)
- [`pcrbias_mock`](https://jsilve24.github.io/fido/reference/pcrbias_mock.md)
  : Data from Silverman et al. (2019) bioRxiv
- [`pibble_sim()`](https://jsilve24.github.io/fido/reference/pibble_sim.md)
  : Simulate simple pibble dataset and priors (for testing)
- [`RISK_CCFA`](https://jsilve24.github.io/fido/reference/RISK_CCFA.md)
  : Data from Gevers et al. (2014)
- [`RISK_CCFA`](https://jsilve24.github.io/fido/reference/RISK_CCFA_otu.md)
  : Data from Gevers et al. (2014)
- [`RISK_CCFA`](https://jsilve24.github.io/fido/reference/RISK_CCFA_sam.md)
  : Data from Gevers et al. (2014)
- [`RISK_CCFA`](https://jsilve24.github.io/fido/reference/RISK_CCFA_tax.md)
  : Data from Gevers et al. (2014)
- [`Y`](https://jsilve24.github.io/fido/reference/Y.md) : Data from
  Silverman et al. (2019) bioRxiv

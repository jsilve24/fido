% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bassetfit_s3.R
\name{predict.bassetfit}
\alias{predict.bassetfit}
\title{Predict using basset}
\usage{
\method{predict}{bassetfit}(
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
}
\arguments{
\item{object}{An object of class pibblefit}

\item{newdata}{An optional matrix for which to evaluate prediction.}

\item{response}{Options = "Lambda":Mean of regression, "Eta", "Y": counts}

\item{size}{the number of counts per sample if response="Y" (as vector or matrix), 
default if newdata=NULL and response="Y" is to use colsums of m$Y. Otherwise
uses median colsums of object$Y as default. If passed as a matrix should have dimensions
ncol(newdata) x iter.}

\item{use_names}{if TRUE apply names to output}

\item{summary}{if TRUE, posterior summary of predictions are returned rather
than samples}

\item{iter}{number of iterations to return if NULL uses object$iter}

\item{from_scratch}{should predictions of Y come from fitted Eta or from 
predictions of Eta from posterior of Lambda? (default: false)}

\item{...}{other arguments passed to summarise_posterior}
}
\value{
(if summary==FALSE) array D x N x iter; (if summary==TRUE) 
tibble with calculated posterior summaries
}
\description{
Predict using basset
}
\details{
currently only implemented for pibblefit objects in coord_system "default"
"alr", or "ilr".
}

\name{correctReadCounts}
\alias{correctReadCounts}
\title{Correct read counts for GC-content bias}
\usage{
correctReadCounts(RCA, cores = 1L, verbose = TRUE)
}
\arguments{
\item{RCA}{An object of class \code{RCA}, as returned by \code{RCA()}.}

\item{cores}{Number of cores used for zero-inflated negative binomial model
fitting via \code{glmmTMB} parallel control (default \code{1L}). Values
outside \code{1} to the number of detected cores fall back to \code{1L}.}

\item{verbose}{Logical; if \code{TRUE} (default), print progress messages.}
}
\value{
The input \code{RCA} object with:
\describe{
  \item{\code{data}}{Updated with logical columns \code{valid} and
    \code{ideal}, and numeric column \code{cor.gc} (GC-corrected coverage).}
  \item{\code{corrected}}{Set to \code{TRUE}.}
  \item{\code{fit}}{The fitted \code{glmmTMB} model object.}
  \item{\code{outliers}}{List of GC and coverage bounds used for filtering:
    \code{gc_lower_bound}, \code{gc_upper_bound}, and
    \code{cov_upper_bound}.}
}
}
\description{
Filters outlier windows and fits a zero-inflated negative binomial model of
coverage as a smooth function of GC-content (with \code{subgenome} as a
fixed covariate) to produce GC-corrected coverage values.
}
\details{
Outlier quantiles are fixed at 1\% / 99\% for GC content and at the 99\%
quantile for coverage.

Correction procedure:
\enumerate{
  \item Mark windows with non-negative coverage and positive GC content as
    \code{valid}.
  \item Among valid windows, mark as \code{ideal} those whose GC content
    lies between the 1\% and 99\% quantiles and whose coverage is at or
    below the 99\% quantile. These windows are used for model fitting.
  \item Fit a zero-inflated negative binomial model on ideal windows, using
    GC-content as a smooth term and subgenome as a fixed covariate. The same
    terms are used for the zero-inflation component.
  \item Predict expected coverage for every window at its observed GC, and
    again at a reference GC equal to the median GC of ideal windows.
  \item Derive GC-corrected coverage by scaling the observed coverage by the
    ratio of reference-expected to fitted coverage.
  \item Further restrict \code{ideal} to windows whose corrected coverage is
    below the 99\% quantile of \code{cor.gc}.
}

Stops with an error if fewer than 2 valid windows remain, if fewer than 10
ideal windows remain before model fitting, or if no ideal windows remain
after the post-correction coverage filter.
}
\examples{
\dontrun{
rca <- correctReadCounts(rca)
rca <- correctReadCounts(rca, cores = 4L)
}
}

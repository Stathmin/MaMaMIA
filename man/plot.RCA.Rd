\name{plot.RCA}
\alias{plot.RCA}
\title{Plot method for RCA objects}
\usage{
\method{plot}{RCA}(x, plot.type = c("orig_cov", "corr_cov"), show_outliers = FALSE, ...)
}
\arguments{
\item{x}{An object of class \code{RCA}.}

\item{plot.type}{Character string selecting the profile to plot. One of
\code{"orig_cov"} (raw coverage) or \code{"corr_cov"} (GC-corrected
coverage).}

\item{show_outliers}{Logical; if \code{FALSE} (default), only windows marked
as \code{ideal} are plotted. If \code{TRUE}, all windows are shown.}

\item{...}{Additional arguments (ignored); present for compatibility with
the generic \code{plot()}.}
}
\value{
A \code{ggplot2::ggplot()} object.
}
\description{
Draws a genome-wide coverage profile for an \code{RCA} object, either as raw
coverage or as GC-corrected coverage.
}
\details{
If GC-bias correction has not been performed (\code{x$corrected} is
\code{FALSE}), two warnings are issued and the method falls back to
\code{plot.type = "orig_cov"} with \code{show_outliers = TRUE}. Stops if no
windows remain after filtering.
}
\examples{
\dontrun{
plot(rca, plot.type = "orig_cov")
plot(rca, plot.type = "corr_cov")
plot(rca, plot.type = "orig_cov", show_outliers = TRUE)
}
}

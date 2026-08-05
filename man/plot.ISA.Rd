\name{plot.ISA}
\alias{plot.ISA}
\title{Plot method for ISA objects}
\usage{
\method{plot}{ISA}(x, plot.type = c("mirror", "diff"), all_segments = FALSE, ...)
}
\arguments{
\item{x}{An object of class \code{ISA}, as returned by \code{segments()}.}

\item{plot.type}{Character string selecting the plot style. One of
\code{"mirror"} (donor coverage above zero and recipient coverage below) or
\code{"diff"} (coverage difference with horizontal segment summary metrics).}

\item{all_segments}{Logical; if \code{FALSE} (default), only segments marked
\code{putative_introgression} are drawn on difference plots. If
\code{TRUE}, all segments are shown.}

\item{...}{Additional arguments (ignored); present for compatibility with
the generic \code{plot()}.}
}
\value{
A \code{ggplot2::ggplot()} object.
}
\description{
Draws pairwise donor-recipient coverage profiles from an \code{ISA} object,
either as a mirrored plot of GC-corrected coverages or as a difference
profile with segment summaries (mean or median, per \code{x$meta$metrics})
overlaid.
}
\examples{
\dontrun{
plot(isa, plot.type = "mirror")
plot(isa, plot.type = "diff")
plot(isa, plot.type = "diff", all_segments = TRUE)
}
}

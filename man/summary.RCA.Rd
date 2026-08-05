\name{summary.RCA}
\alias{summary.RCA}
\alias{print.summary.RCA}
\title{Summary method for RCA objects}
\usage{
\method{summary}{RCA}(object, ...)

\method{print}{summary.RCA}(x, ...)
}
\arguments{
\item{object}{An object of class \code{RCA}.}

\item{x}{An object of class \code{summary.RCA}.}

\item{...}{Additional arguments (ignored); present for compatibility with
the generics.}
}
\value{
An object of class \code{summary.RCA}: a list with
\code{n_chromosomes}, \code{n_subgenomes},
\code{chromosomes_per_subgenome} (named integer vector of chromosome counts
per subgenome label), \code{corrected}, and \code{outliers} (the GC and
coverage bounds from \code{correctReadCounts()}, or \code{NULL} if
correction has not been applied). Printing the object shows a
human-readable summary, including outlier thresholds when available.
}
\description{
Builds a structured summary of chromosome and subgenome counts for an
\code{RCA} object, and whether GC-bias correction has been applied.
}
\examples{
\dontrun{
summary(rca)
s <- summary(rca)
s$n_chromosomes
s$outliers
}
}

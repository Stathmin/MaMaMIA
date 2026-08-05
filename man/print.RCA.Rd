\name{print.RCA}
\alias{print.RCA}
\title{Print method for RCA objects}
\usage{
\method{print}{RCA}(x, ...)
}
\arguments{
\item{x}{An object of class \code{RCA}.}

\item{...}{Additional arguments (ignored); present for compatibility with
the generic \code{print()}.}
}
\value{
\code{x}, invisibly.
}
\description{
Prints a one-line overview of an \code{RCA} object: total chromosome count,
donor and recipient counts, window count, and whether GC-bias correction
has been applied.
}
\examples{
\dontrun{
print(rca)
}
}

\name{print.ISA}
\alias{print.ISA}
\title{Print method for ISA objects}
\usage{
\method{print}{ISA}(x, ...)
}
\arguments{
\item{x}{An object of class \code{ISA}, as returned by \code{segments()}.}

\item{...}{Additional arguments (ignored); present for compatibility with
the generic \code{print()}.}
}
\value{
\code{x}, invisibly.
}
\description{
Prints a one-line overview of an \code{ISA} object (pair count, segment
count, putative introgression count, and active metric).
}
\examples{
\dontrun{
print(isa)
}
}

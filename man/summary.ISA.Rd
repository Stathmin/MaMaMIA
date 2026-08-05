\name{summary.ISA}
\alias{summary.ISA}
\alias{print.summary.ISA}
\title{Summary method for ISA objects}
\usage{
\method{summary}{ISA}(object, unit = c("bp", "Kb", "Mb", "Gb"), digits = 2, ...)

\method{print}{summary.ISA}(x, ...)
}
\arguments{
\item{object}{An object of class \code{ISA}, as returned by \code{segments()}.}

\item{unit}{Length unit for reported coordinates and fragment lengths. One of
\code{"bp"} (default), \code{"Kb"}, \code{"Mb"}, or \code{"Gb"}.}

\item{digits}{Non-negative integer giving the number of decimal places used
when scaling coordinates and lengths to \code{unit} (default \code{2}).}

\item{x}{An object of class \code{summary.ISA}.}

\item{...}{Additional arguments (ignored); present for compatibility with
the generics.}
}
\value{
An object of class \code{summary.ISA}: a list with
\code{introgressions} (a data frame of putative introgressions with donor
and recipient coordinates and lengths), \code{n_introgressions},
\code{unit}, and \code{digits}. Printing the object formats coordinates
in the requested unit.
}
\description{
Builds a structured summary of putative introgressions from an \code{ISA}
object. Coordinates and fragment lengths are taken from the donor and
recipient windows spanning each flagged segment.
}
\details{
Only segments with \code{putative_introgression = TRUE} are included. When
no such segments exist, printing reports that none were found. Each
introgression is shown as a two-row table (donor and recipient) with
chromosome, subgenome, From, To, and Length columns.
}
\examples{
data(triticum)
rca <- RCA(don_cov, don_gc, rec_cov, rec_gc, meta)
\donttest{
## Note: correctReadCounts() can take a while on the full dataset
rca <- correctReadCounts(rca, cores = 4L, verbose = FALSE)
}

## For demonstration purposes, assume that GC bias correction was performed, 
## by copying the 'cov' column to 'cor.gc'.
rca$data$valid <- TRUE
gc_q <- quantile(rca$data$gc[rca$data$valid], c(0.01, 0.99))
cov_q <- quantile(rca$data$cov[rca$data$valid], 0.99)
rca$data$ideal <- rca$data$valid & 
                    (rca$data$gc >= gc_q[1]) & 
                    (rca$data$gc <= gc_q[2]) & 
                    (rca$data$cov <= cov_q)
rca$corrected <- TRUE
rca$outliers <- list(
        gc_lower_bound = gc_q[1],
        gc_upper_bound = gc_q[2],
        cov_upper_bound = cov_q
    )
rca$data$cor.gc <- rca$data$cov

isa <- segments(
  rca,
  target_pairs = c("OY997261.1" = "CP155614.1"),
  verbose = FALSE
)
summary(isa)
s <- summary(isa, unit = "Mb", digits = 3)
s$n_introgressions
}

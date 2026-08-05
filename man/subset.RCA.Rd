\name{subset.RCA}
\alias{subset.RCA}
\title{Subset method for RCA objects}
\usage{
\method{subset}{RCA}(
  x,
  don_chromlist = NULL,
  rec_chromlist = NULL,
  don_subgenomes = NULL,
  rec_subgenomes = NULL,
  ...
)
}
\arguments{
\item{x}{An object of class \code{RCA}.}

\item{don_chromlist}{Character or numeric vector of donor chromosome IDs to
retain, or \code{NULL} (default) to leave the donor side unchanged unless
\code{don_subgenomes} is set. Empty vectors are not allowed.}

\item{rec_chromlist}{Character or numeric vector of recipient chromosome IDs
to retain, or \code{NULL} (default) to leave the recipient side unchanged
unless \code{rec_subgenomes} is set. Empty vectors are not allowed.}

\item{don_subgenomes}{Character vector of donor subgenome labels whose
chromosomes should be retained, or \code{NULL} (default). When provided
together with \code{don_chromlist}, matching chromosomes are combined
(union). Labels must already appear in \code{x$meta$meta$subgenome}.}

\item{rec_subgenomes}{Character vector of recipient subgenome labels whose
chromosomes should be retained, or \code{NULL} (default). Combined with
\code{rec_chromlist} as for the donor side.}

\item{...}{Additional arguments (ignored); present for compatibility with
the generic \code{subset()}.}
}
\value{
An \code{RCA} object containing only the selected chromosomes, with
updated \code{data}, \code{meta$meta}, \code{meta$don_chr_ids}, and
\code{meta$rec_chr_ids}. Other components (\code{corrected}, \code{fit},
\code{outliers}) are left unchanged.
}
\description{
Restricts an \code{RCA} object to selected donor and/or recipient chromosomes,
specified either by chromosome ID lists or by subgenome labels. At least one
of the chromosome or subgenome arguments must be non-\code{NULL}. Omitting
filters for one side leaves that side unchanged. The subset must leave at
least one donor and one recipient chromosome.
}
\examples{
data(triticum)
rca <- RCA(don_cov, don_gc, rec_cov, rec_gc, meta)
rca_sub <- subset(rca, don_subgenomes = "At", rec_subgenomes = "A")
rca_sub <- subset(rca, don_subgenomes = "G")
rca_sub <- subset(
  rca,
  don_chromlist = "OY997261.1",
  rec_chromlist = "CP155614.1"
)
rca_sub
}

\docType{data}
\name{triticum}
\alias{triticum}
\alias{don_cov}
\alias{don_gc}
\alias{rec_cov}
\alias{rec_gc}
\alias{meta}
\title{Triticum introgression example dataset}
\format{
A dataset with five data frames:

\describe{
  \item{\code{don_cov}}{Donor coverage (\emph{Triticum timopheevii}):
    a data frame with 14000 rows and 4 columns (no column names;
    columns are positional):
    \describe{
      \item{V1}{Chromosome accession ID (14 chromosomes; At and G subgenomes).}
      \item{V2}{Interval start coordinate (0-based).}
      \item{V3}{Interval end coordinate (exclusive).}
      \item{V4}{Number of reads mapped in the interval.}
    }}
  \item{\code{don_gc}}{Donor GC-content (\emph{Triticum timopheevii}):
    a data frame with 14000 rows and 4 columns, in the same windows as
    \code{don_cov}:
    \describe{
      \item{V1}{Chromosome accession ID (match \code{don_cov}).}
      \item{V2}{Interval start coordinate (0-based).}
      \item{V3}{Interval end coordinate (exclusive).}
      \item{V4}{GC-content of the interval (proportion between 0 and 1).}
    }}
  \item{\code{rec_cov}}{Recipient coverage (\emph{Triticum aestivum}
    Chinese Spring T2T): a data frame with 21000 rows and 4 columns
    (no column names; columns are positional):
    \describe{
      \item{V1}{Chromosome accession ID (21 chromosomes; A, B, and D
        subgenomes).}
      \item{V2}{Interval start coordinate (0-based).}
      \item{V3}{Interval end coordinate (exclusive).}
      \item{V4}{Number of reads mapped in the interval.}
    }}
  \item{\code{rec_gc}}{Recipient GC-content (\emph{Triticum aestivum}
    Chinese Spring T2T): a data frame with 21000 rows and 4 columns,
    in the same windows as \code{rec_cov}:
    \describe{
      \item{V1}{Chromosome accession ID (match \code{rec_cov}).}
      \item{V2}{Interval start coordinate (0-based).}
      \item{V3}{Interval end coordinate (exclusive).}
      \item{V4}{GC-content of the interval (proportion between 0 and 1).}
    }}
  \item{\code{meta}}{Chromosome metadata for donor and recipient genomes:
    a data frame with 35 rows and 3 columns:
    \describe{
      \item{V1}{Chromosome accession ID (matches IDs in coverage and
        GC-content tables).}
      \item{V2}{Short chromosome name (e.g. \code{"1A"}, \code{"1At"}, \code{"1G"}).}
      \item{V3}{Subgenome label: \code{"A"}, \code{"B"}, \code{"D"}, \code{"At"}, or \code{"G"}
        (7 chromosomes each).}
    }}
}
}
\source{
Internally generated.
}
\usage{
data(triticum)
}
\description{
Example windowed coverage, GC-content, and chromosome metadata for
introgression analysis between donor \emph{Triticum timopheevii} and
recipient \emph{Triticum aestivum} Chinese Spring T2T. Intended for use
with \code{RCA()} as \code{dcov}, \code{dgc}, \code{rcov}, \code{rgc}, and
\code{meta}.
}
\keyword{datasets}

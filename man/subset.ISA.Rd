\name{subset.ISA}
\alias{subset.ISA}
\title{Subset method for ISA objects}
\usage{
\method{subset}{ISA}(
  x,
  min_width = NULL,
  min_segment_median = NULL,
  min_segment_mean = NULL,
  threshold_type = NULL,
  target_pairs = NULL,
  ...
)
}
\arguments{
\item{x}{An object of class \code{ISA}, as returned by \code{segments()}.}

\item{min_width}{New minimum segment width in windows, or \code{NULL}
(default) to leave unchanged. Must be \code{>= 2}. When supplied, segments
shorter than this width are dropped and \code{x$param$min.width} is
updated.}

\item{min_segment_median}{New minimum segment median threshold, or
\code{NULL} (default) to leave unchanged. Requires that the ISA metrics
column is \code{seg.median}. Cannot be set together with
\code{min_segment_mean}.}

\item{min_segment_mean}{New minimum segment mean threshold, or
\code{NULL} (default) to leave unchanged. Requires that the ISA metrics
column is \code{seg.mean}. Cannot be set together with
\code{min_segment_median}.}

\item{threshold_type}{Optional character string \code{"mean"} or
\code{"median"} to switch which segment summary statistic is stored in
\code{x$meta$metrics}. Switching to a different metric than the current
one requires the matching \code{min_segment_mean} or
\code{min_segment_median} argument so the threshold is not reused across
incomparable statistics.}

\item{target_pairs}{Optional named vector mapping donor chromosome IDs
(\code{names}) to recipient chromosome IDs (\code{values}). When supplied,
restricts \code{data}, \code{out}, \code{segments}, and metadata to those
pairs. \code{don_chr_ids} / \code{rec_chr_ids} are set to the unique
chromosome IDs involved. All IDs must already be present in the ISA object.}

\item{...}{Additional arguments (ignored); present for compatibility with
the generic \code{subset()}.}
}
\value{
An \code{ISA} object with updated threshold metadata, minimum width,
and/or restricted chromosome pairs. The \code{putative_introgression}
column is recalculated from the active metrics column and threshold.
Window data, pairwise output, and segment tables are filtered when
\code{target_pairs} is supplied. Segments shorter than \code{min_width}
are dropped when that argument is set.
}
\description{
Updates threshold settings, minimum segment width, and/or restricts an
\code{ISA} object to selected donor-recipient chromosome pairs. Does not
re-run segmentation.
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
  target_pairs = c(
    "OY997261.1" = "CP155614.1",
    "OY997263.1" = "CP155611.1"
  ),
  verbose = FALSE
)
isa_sub <- subset(isa, target_pairs = c("OY997261.1" = "CP155614.1"))
isa_sub <- subset(isa, threshold_type = "mean", min_segment_mean = 5)
isa_sub <- subset(isa, min_width = 10)
isa_sub
}

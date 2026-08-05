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
data(triticum)
rca <- RCA(don_cov, don_gc, rec_cov, rec_gc, meta)
plot(rca, plot.type = "orig_cov")
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

plot(rca, plot.type = "corr_cov")
plot(rca, plot.type = "orig_cov", show_outliers = TRUE)
}

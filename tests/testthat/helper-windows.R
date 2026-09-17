# Shared fixtures for the test suite.
#
# Everything here is deliberately tiny and free of model fitting so the default
# test run stays fast enough for CRAN. The expensive glmmTMB paths are exercised
# separately in test-correction.R under skip_on_cran().

# Contiguous half-open windows covering one chromosome.
make_windows <- function(chr, n, value, width = 100) {
    start <- seq(0, by = width, length.out = n)
    data.frame(
        chr_id = chr,
        start = start,
        stop = start + width,
        value = value
    )
}

make_meta <- function(chr_id, chr_name, subgenome) {
    data.frame(chr_id = chr_id, chr_name = chr_name, subgenome = subgenome)
}

# Windowed tables for two donor and two recipient chromosomes.
small_inputs <- function(n = 40L) {
    cov <- rep(10, n)
    gc <- seq(0.3, 0.6, length.out = n)
    list(
        dcov = rbind(make_windows("d1", n, cov), make_windows("d2", n, cov)),
        dgc = rbind(make_windows("d1", n, gc), make_windows("d2", n, gc)),
        rcov = rbind(make_windows("r1", n, cov), make_windows("r2", n, cov)),
        rgc = rbind(make_windows("r1", n, gc), make_windows("r2", n, gc)),
        meta = make_meta(
            c("d1", "d2", "r1", "r2"),
            c("1At", "2At", "1A", "2A"),
            c("At", "G", "A", "B")
        )
    )
}

small_rca <- function(n = 40L) {
    do.call(RCA, small_inputs(n))
}

# A single donor/recipient pair in which the donor gains coverage over the second
# half of the chromosome, so CBS has one obvious change point to find. The donor
# profile is deliberately not constant: it is asymmetric, so reversing the
# chromosome visibly changes the difference profile.
step_inputs <- function(n = 60L) {
    donor_cov <- rep(8, n)
    donor_cov[31:n] <- 16
    recipient_cov <- rep(10, n)
    gc <- seq(0.35, 0.55, length.out = n)
    list(
        dcov = make_windows("d1", n, donor_cov),
        dgc = make_windows("d1", n, gc),
        rcov = make_windows("r1", n, recipient_cov),
        rgc = make_windows("r1", n, gc),
        meta = make_meta(c("d1", "r1"), c("1At", "1A"), c("At", "A"))
    )
}

# An RCA object flagged as GC-corrected without fitting a model, mirroring the
# shortcut used in the package vignette's segmentation examples.
step_rca <- function(n = 60L) {
    x <- do.call(RCA, step_inputs(n))
    x$data$valid <- TRUE
    x$data$ideal <- TRUE
    x$data$cor.gc <- x$data$cov
    x$corrected <- TRUE
    x$outliers <- list(
        gc_lower_bound = 0.01,
        gc_upper_bound = 0.99,
        cov_upper_bound = Inf
    )
    x
}

# One donor and one recipient chromosome taken from the packaged example data,
# used by the tests that need a well-posed model fit. The packaged tables carry
# positional column names (V1, V2, ...), so they are selected by position.
triticum_one_pair <- function() {
    utils::data("triticum", package = "MaMaMIA", envir = environment())

    pick <- function(df, id, value_col) {
        out <- df[df[[1L]] == id, , drop = FALSE]
        names(out) <- c("chr_id", "start", "stop", value_col)
        out
    }

    don_id <- unique(don_cov[[1L]])[1L]
    rec_id <- unique(rec_cov[[1L]])[1L]
    keep <- c(don_id, rec_id)

    list(
        dcov = pick(don_cov, don_id, "cov"),
        dgc = pick(don_gc, don_id, "gc"),
        rcov = pick(rec_cov, rec_id, "cov"),
        rgc = pick(rec_gc, rec_id, "gc"),
        meta = stats::setNames(
            meta[meta[[1L]] %in% keep, , drop = FALSE],
            c("chr_id", "chr_name", "subgenome")
        )
    )
}

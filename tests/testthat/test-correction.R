# These tests fit the zero-inflated negative binomial model, which takes tens of
# seconds, so they are skipped on CRAN. The cheap validation of the surrounding
# pipeline lives in the other test files.

test_that("correctReadCounts() stores GC-corrected coverage", {
    skip_on_cran()

    x <- do.call(RCA, triticum_one_pair())
    res <- correctReadCounts(x, cores = 1L, verbose = FALSE)

    expect_s3_class(res, "RCA")
    expect_true(res$corrected)
    expect_s3_class(res$fit, "glmmTMB")
    expect_true(all(c("valid", "ideal", "cor.gc") %in% names(res$data)))
    expect_type(res$data$valid, "logical")
    expect_type(res$data$ideal, "logical")
    expect_true(is.numeric(res$data$cor.gc))
    expect_true(all(is.finite(res$data$cor.gc)))
    expect_true(any(res$data$ideal))

    expect_length(res$outliers, 3L)
    expect_named(
        res$outliers,
        c("gc_lower_bound", "gc_upper_bound", "cov_upper_bound")
    )
    expect_true(res$outliers$gc_lower_bound < res$outliers$gc_upper_bound)
})

test_that("correctReadCounts() reports the correction in summary()", {
    skip_on_cran()

    x <- do.call(RCA, triticum_one_pair())
    res <- correctReadCounts(x, cores = 1L, verbose = FALSE)
    s <- summary(res)

    expect_true(s$corrected)
    expect_false(is.null(s$outliers))
    expect_output(print(s), "GC-bias corrected: Yes")
})

test_that("correctReadCounts() validates its input", {
    expect_error(correctReadCounts(list()), "not RCA")
})

test_that("a cached fit is reused and returns identical results", {
    skip_on_cran()

    x <- do.call(RCA, triticum_one_pair())
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path), add = TRUE)

    cold <- correctReadCounts(x, cores = 1L, verbose = FALSE, cache = path)
    expect_true(file.exists(path))

    warm <- correctReadCounts(x, cores = 1L, verbose = FALSE, cache = path)
    expect_identical(warm$data$cor.gc, cold$data$cor.gc)
    expect_identical(warm$data$ideal, cold$data$ideal)
    expect_identical(glmmTMB::fixef(warm$fit), glmmTMB::fixef(cold$fit))

    # A stale entry is not reused: the fitting data changed.
    reversed <- reverseWindows(x, chr_ids = unique(x$data$chr_id)[1L])
    refitted <- correctReadCounts(reversed, cores = 1L, verbose = FALSE, cache = path)
    expect_true(refitted$corrected)
})

test_that("the fast ZINB prediction reproduces glmmTMB::predict()", {
    skip_on_cran()

    x <- do.call(RCA, triticum_one_pair())
    res <- correctReadCounts(x, cores = 1L, verbose = FALSE)
    gc_ref <- stats::median(res$data$gc[res$data$ideal], na.rm = TRUE)

    fast <- MaMaMIA:::predict_zinb_response(
        res$fit,
        gc = res$data$gc,
        subgenome = res$data$subgenome,
        gc_ref = gc_ref
    )

    nd <- data.frame(gc = res$data$gc, subgenome = res$data$subgenome)
    reference_actual <- stats::predict(res$fit, newdata = nd, type = "response")
    reference_ref <- stats::predict(
        res$fit,
        newdata = transform(nd, gc = gc_ref),
        type = "response"
    )

    expect_equal(fast$actual, as.numeric(reference_actual), tolerance = 1e-8)
    expect_equal(fast$ref, as.numeric(reference_ref), tolerance = 1e-8)
    expect_length(fast$actual, nrow(res$data))
    expect_length(fast$ref, nrow(res$data))

    # Prediction at the reference GC depends only on the subgenome.
    spread <- tapply(fast$ref, res$data$subgenome, function(z) diff(range(z)))
    expect_true(all(unlist(spread) < 1e-8))
})

test_that("the fast ZINB prediction falls back to predict() when unavailable", {
    skip_on_cran()

    x <- do.call(RCA, triticum_one_pair())
    res <- correctReadCounts(x, cores = 1L, verbose = FALSE)
    gc_ref <- stats::median(res$data$gc[res$data$ideal], na.rm = TRUE)

    # Removing the parameter vector the fast path relies on must not change the
    # values returned, it should just take the slower route.
    broken <- res$fit
    broken$fit$parfull <- NULL

    fast <- MaMaMIA:::predict_zinb_response(
        res$fit,
        gc = res$data$gc,
        subgenome = res$data$subgenome,
        gc_ref = gc_ref
    )
    fallback <- MaMaMIA:::predict_zinb_response(
        broken,
        gc = res$data$gc,
        subgenome = res$data$subgenome,
        gc_ref = gc_ref
    )

    expect_equal(fast$actual, fallback$actual, tolerance = 1e-8)
    expect_equal(fast$ref, fallback$ref, tolerance = 1e-8)
})

test_that("correctReadCounts() output feeds segmentation", {
    skip_on_cran()

    x <- do.call(RCA, triticum_one_pair())
    res <- correctReadCounts(x, cores = 1L, verbose = FALSE)
    pairs <- stats::setNames(res$meta$rec_chr_ids, res$meta$don_chr_ids)

    isa <- segments(res, target_pairs = pairs, verbose = FALSE)

    expect_s3_class(isa, "ISA")
    expect_true(all(isa$segments$chr_id.don %in% res$meta$don_chr_ids))
    expect_true(all(isa$segments$chr_id.rec %in% res$meta$rec_chr_ids))
    expect_true(all(is.finite(isa$segments$seg.median)))
})

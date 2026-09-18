test_that("the fit cache key follows the windows, cores and versions", {
    windows <- data.frame(
        chr_id = "c1",
        iid = 1:4,
        cov = c(10, 11, 12, 13),
        gc = c(0.40, 0.41, 0.42, 0.43),
        subgenome = "A"
    )
    call <- MaMaMIA:::coverage_model_call(4L)

    key <- MaMaMIA:::fit_cache_key(windows, 4L, call)
    expect_type(key, "character")
    expect_length(key, 1L)
    expect_identical(key, MaMaMIA:::fit_cache_key(windows, 4L, call))

    # Anything that changes the fitted values must change the key.
    expect_false(identical(
        key,
        MaMaMIA:::fit_cache_key(windows, 8L, MaMaMIA:::coverage_model_call(8L))
    ))
    expect_false(identical(
        key,
        MaMaMIA:::fit_cache_key(transform(windows, cov = cov + 1), 4L, call)
    ))
    expect_false(identical(
        key,
        MaMaMIA:::fit_cache_key(transform(windows, gc = gc + 1e-3), 4L, call)
    ))
    expect_false(identical(
        key,
        MaMaMIA:::fit_cache_key(transform(windows, subgenome = "B"), 4L, call)
    ))
    expect_false(identical(
        key,
        MaMaMIA:::fit_cache_key(windows[-1L, ], 4L, call)
    ))
})

test_that("the fit cache key follows the model specification", {
    windows <- data.frame(
        chr_id = "c1",
        iid = 1:4,
        cov = c(10, 11, 12, 13),
        gc = c(0.40, 0.41, 0.42, 0.43),
        subgenome = "A"
    )
    call <- MaMaMIA:::coverage_model_call(4L)
    key <- MaMaMIA:::fit_cache_key(windows, 4L, call)

    # A rewritten model must never be served a fit from the old one: the package
    # version alone cannot catch this, since it does not move during development.
    wider <- as.list(call)
    wider[[2L]] <- quote(cov ~ s(gc, k = 15) + subgenome)
    expect_false(identical(
        key,
        MaMaMIA:::fit_cache_key(windows, 4L, as.call(wider))
    ))

    no_zi <- as.list(call)
    no_zi$ziformula <- quote(~ 1)
    expect_false(identical(
        key,
        MaMaMIA:::fit_cache_key(windows, 4L, as.call(no_zi))
    ))

    # The data argument itself is not part of the specification: the windows are
    # keyed separately, and by the columns the model uses rather than all of them.
    relabelled <- as.list(call)
    relabelled$data <- quote(anything)
    expect_identical(
        key,
        MaMaMIA:::fit_cache_key(windows, 4L, as.call(relabelled))
    )
})

test_that("cached fits are reused only when the key matches", {
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path), add = TRUE)

    # No cache argument, or no file yet.
    expect_null(MaMaMIA:::read_cached_fit(NULL, "key"))
    expect_null(MaMaMIA:::read_cached_fit(path, "key"))

    fake <- structure(list(fixef = 1), class = "glmmTMB")
    MaMaMIA:::write_cached_fit(path, "key", fake)
    expect_true(file.exists(path))
    expect_identical(MaMaMIA:::read_cached_fit(path, "key"), fake)
    expect_identical(
        readRDS(path)$format,
        MaMaMIA:::CACHE_LAYOUT_VERSION
    )

    # A different key must not reuse the entry.
    expect_null(MaMaMIA:::read_cached_fit(path, "other"))

    # A corrupt cache is ignored rather than throwing.
    writeLines("not an rds file", path)
    expect_null(MaMaMIA:::read_cached_fit(path, "key"))

    # So is one written by a different layout, even with a matching key.
    for (stale in list(NULL, MaMaMIA:::CACHE_LAYOUT_VERSION + 1L)) {
        entry <- list(format = stale, key = "key", fit = fake)
        saveRDS(entry, path)
        expect_null(MaMaMIA:::read_cached_fit(path, "key"))
    }

    # With caching off nothing is written.
    MaMaMIA:::write_cached_fit(NULL, "key", fake)
    expect_false(file.exists(paste0(path, ".tmp")))
})

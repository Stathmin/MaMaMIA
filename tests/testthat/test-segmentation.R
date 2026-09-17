test_that("segments() returns an ISA holding the requested pair", {
    x <- step_rca()
    isa <- segments(x, target_pairs = c("d1" = "r1"), verbose = FALSE)

    expect_s3_class(isa, "ISA")
    expect_true(nrow(isa$segments) >= 1L)
    expect_true(all(c(
        "chr_id.don", "chr_id.rec", "loc.start", "loc.end", "num.mark",
        "seg.mean", "seg.sd", "seg.median", "seg.mad",
        "putative_introgression",
        "chr_name.don", "subgenome.don", "chr_name.rec", "subgenome.rec"
    ) %in% names(isa$segments)))

    expect_type(isa$segments$putative_introgression, "logical")
    expect_true(all(isa$segments$chr_id.don == "d1"))
    expect_true(all(isa$segments$chr_id.rec == "r1"))
    expect_true(all(isa$segments$loc.start >= 1))
    expect_true(all(isa$segments$loc.end <= 60))

    # The recipient loses coverage over the second half of the chromosome, so at
    # least one segment must show a donor excess.
    expect_true(any(isa$segments$seg.median > 0))

    expect_equal(isa$meta$metrics, "seg.median")
    expect_setequal(isa$meta$don_chr_ids, "d1")
    expect_setequal(isa$meta$rec_chr_ids, "r1")
})

test_that("segments() honours the introgression threshold", {
    x <- step_rca()

    permissive <- segments(x, target_pairs = c("d1" = "r1"), verbose = FALSE)
    strict <- segments(
        x,
        target_pairs = c("d1" = "r1"),
        min_segment_median = 1e9,
        verbose = FALSE
    )

    expect_true(any(permissive$segments$putative_introgression))
    expect_false(any(strict$segments$putative_introgression))
    expect_equal(nrow(strict$segments), nrow(permissive$segments))
})

test_that("segments() drops segments shorter than min_width", {
    x <- step_rca()

    isa <- segments(
        x,
        target_pairs = c("d1" = "r1"),
        min_width = 1000L,
        verbose = FALSE
    )

    expect_equal(nrow(isa$segments), 0L)
})

test_that("segments() validates its inputs", {
    corrected <- step_rca()

    expect_error(segments(1), "not RCA")
    expect_error(
        segments(do.call(RCA, step_inputs()), target_pairs = c("d1" = "r1")),
        "corrected"
    )
    expect_error(
        segments(corrected, target_pairs = c("nope" = "r1"), verbose = FALSE),
        "not present"
    )
    expect_error(
        segments(corrected, target_pairs = c("d1" = "nope"), verbose = FALSE),
        "not present"
    )
})

test_that("segments() is reproducible regardless of the caller's RNG state", {
    x <- step_rca()
    tp <- c("d1" = "r1")

    set.seed(999)
    first <- segments(x, target_pairs = tp, verbose = FALSE)
    set.seed(123)
    second <- segments(x, target_pairs = tp, verbose = FALSE)

    expect_identical(first$segments, second$segments)
    expect_identical(
        first$segments,
        segments(x, target_pairs = tp, verbose = FALSE)$segments
    )
})

test_that("segments() leaves the caller's RNG stream untouched", {
    x <- step_rca()

    set.seed(7)
    before <- stats::runif(3)
    set.seed(7)
    invisible(segments(x, target_pairs = c("d1" = "r1"), verbose = FALSE))

    expect_identical(before, stats::runif(3))
})

test_that("segments() stores and validates its seed", {
    x <- step_rca()

    expect_equal(
        segments(x, target_pairs = c("d1" = "r1"), seed = 5L, verbose = FALSE)$param$seed,
        5L
    )
    expect_equal(
        segments(x, target_pairs = c("d1" = "r1"), verbose = FALSE)$param$seed,
        1L
    )
    expect_error(
        segments(x, target_pairs = c("d1" = "r1"), seed = "x", verbose = FALSE)
    )
    expect_error(
        segments(x, target_pairs = c("d1" = "r1"), seed = c(1L, 2L), verbose = FALSE)
    )
})

test_that("segmentation parameters are validated", {
    x <- step_rca()
    tp <- c("d1" = "r1")

    expect_error(segments(x, target_pairs = tp, alpha = 0, verbose = FALSE))
    expect_error(segments(x, target_pairs = tp, alpha = 2, verbose = FALSE))
    expect_error(segments(x, target_pairs = tp, min_width = 1, verbose = FALSE))
    expect_error(segments(x, target_pairs = tp, undo_SD = 0, verbose = FALSE))
    expect_error(segments(x, target_pairs = tp, undo_SD = 11, verbose = FALSE))
})

test_that("reverseWindows() re-segments affected pairs of an ISA", {
    x <- step_rca()
    isa <- segments(x, target_pairs = c("d1" = "r1"), verbose = FALSE)

    reversed <- reverseWindows(isa, chr_ids = "d1")

    expect_s3_class(reversed, "ISA")
    expect_equal(reversed$param, isa$param)
    expect_setequal(reversed$meta$don_chr_ids, "d1")
    expect_setequal(reversed$meta$rec_chr_ids, "r1")
    expect_true(nrow(reversed$segments) >= 1L)

    # The donor profile is asymmetric, so reversing it must change the difference
    # profile that the affected pair is re-segmented on.
    expect_false(identical(reversed$out$diff, isa$out$diff))

    expect_error(reverseWindows(isa, chr_ids = "nope"), "not present")
})

test_that("summary() and print() describe an ISA", {
    isa <- segments(step_rca(), target_pairs = c("d1" = "r1"), verbose = FALSE)

    expect_output(print(isa), "ISA object")
    expect_output(
        print(summary(isa, unit = "Mb")),
        "Putative introgressions found"
    )
})

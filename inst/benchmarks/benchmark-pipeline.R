#!/usr/bin/env Rscript

# Benchmark and consistency harness for the MaMaMIA pipeline.
#
# Times and memory-profiles every step of the documented workflow on the packaged
# `triticum` example, with a piecewise breakdown of correctReadCounts(), which
# dominates the runtime. It also checks that the internals still agree:
#
#   1. the fast ZINB prediction reproduces glmmTMB::predict() on the same fit
#   2. correctReadCounts() is repeatable and its internals sum to the whole call
#   3. segments() is reproducible and leaves the caller's RNG stream untouched
#   4. writeToBedpe() is byte-stable
#   5. with --repeat=N the whole pipeline produces identical fingerprints
#
# Outputs in --out (default "bench-out"): RESULTS.csv (one row per step) and
# FINGERPRINT.txt (environment, output hashes, check results). --compare=<old
# FINGERPRINT.txt> asserts that a change altered no output.
#
# Usage:
#   Rscript inst/benchmarks/benchmark-pipeline.R [--cores=8] [--repeat=2]
#   Rscript inst/benchmarks/benchmark-pipeline.R --compare=bench-out/FINGERPRINT.txt
#
# Exits 1 if a check fails, so it can gate CI.
#
# NOTE: glmmTMB reduces the likelihood in a thread-count-dependent order, so
# fitted values shift by ~1e-7 relative with --cores, and fingerprints are only
# comparable at the same --cores.

# ---------------------------------------------------------------- utilities --

abort <- function(...) stop(..., call. = FALSE)

parse_args <- function(args) {
    opts <- list(cores = 4L, reps = 1L, out = "bench-out", compare = NA_character_)
    for (a in args) {
        if (grepl("^--cores=", a)) {
            opts$cores <- as.integer(sub("^--cores=", "", a))
        } else if (grepl("^--repeat=", a)) {
            opts$reps <- as.integer(sub("^--repeat=", "", a))
        } else if (grepl("^--out=", a)) {
            opts$out <- sub("^--out=", "", a)
        } else if (grepl("^--compare=", a)) {
            opts$compare <- sub("^--compare=", "", a)
        } else {
            abort("unknown argument: ", a, "\nSee the header of this script for usage.")
        }
    }
    if (is.na(opts$cores) || opts$cores < 1L) abort("--cores must be a positive integer")
    if (is.na(opts$reps) || opts$reps < 1L) abort("--repeat must be a positive integer")
    opts
}

# Resident / high-water memory of this process, in kB. Linux only; NA elsewhere.
read_vm <- function(field) {
    status <- "/proc/self/status"
    if (!file.exists(status)) {
        return(NA_real_)
    }
    lines <- readLines(status, warn = FALSE)
    hit <- grep(paste0("^", field, ":"), lines, value = TRUE)
    if (!length(hit)) {
        return(NA_real_)
    }
    as.numeric(sub(" kB", "", sub(paste0("^", field, ":\\s+"), "", hit)))
}

# Peak R heap during the last measurement window, in MB. Ncells are 56 bytes and
# Vcells 8 bytes each on a 64-bit build.
heap_peak_mb <- function(g) {
    g["Ncells", "max used"] * 56 / 1e6 + g["Vcells", "max used"] * 8 / 1e6
}

# Stable content hash, used to compare outputs across runs and code versions.
hash_object <- function(x) {
    f <- tempfile()
    on.exit(unlink(f), add = TRUE)
    saveRDS(x, f, version = 2)
    unname(tools::md5sum(f))
}

new_recorder <- function() {
    e <- new.env(parent = emptyenv())
    e$rows <- list()
    e$checks <- list()
    e
}

record <- function(rec, section, step, elapsed_s, peak_MB, peak_rss_MB, vmhwm_MB, obj_MB) {
    rec$rows[[length(rec$rows) + 1L]] <- data.frame(
        section = section,
        step = step,
        elapsed_s = round(elapsed_s, 3),
        peak_heap_MB = round(peak_MB, 2),
        peak_rss_MB = round(peak_rss_MB, 2),
        vmhwm_MB = round(vmhwm_MB, 2),
        obj_MB = round(obj_MB, 3)
    )
}

# Run `expr`, recording time and memory, and return its value.
#
# Memory is reported three ways, because each alone is misleading: peak_heap_MB
# is the peak R heap during the step (gc() max used), peak_rss_MB the growth of
# the kernel resident high-water mark, and vmhwm_MB that mark once the step is
# done. The growth is 0 whenever the step stayed below an earlier peak, so only
# the absolute column is comparable across steps; and compiled code such as the
# glmmTMB fit allocates outside R's heap, which is why that peak alone
# understates it.
bench <- function(rec, section, step, expr, size = TRUE) {
    hwm_before <- read_vm("VmHWM")
    gc(reset = TRUE)
    timing <- system.time(value <- force(expr))
    g <- gc()
    hwm_after <- read_vm("VmHWM")
    peak_rss <- if (is.na(hwm_before) || is.na(hwm_after)) {
        NA_real_
    } else {
        (hwm_after - hwm_before) / 1024
    }
    obj_mb <- if (size) as.numeric(utils::object.size(value)) / 1e6 else NA_real_
    record(
        rec, section, step, timing[["elapsed"]],
        heap_peak_mb(g), peak_rss, hwm_after / 1024, obj_mb
    )
    invisible(value)
}

check <- function(rec, label, ok, detail = "") {
    rec$checks[[length(rec$checks) + 1L]] <- data.frame(
        check = label,
        result = if (isTRUE(ok)) "PASS" else "FAIL",
        detail = detail
    )
    invisible(ok)
}

print_section <- function(rec, section) {
    rows <- do.call(rbind, rec$rows)
    keep <- rows[rows$section == section, , drop = FALSE]
    if (nrow(keep)) {
        print(keep[, setdiff(names(keep), "section"), drop = FALSE], row.names = FALSE)
    }
}

# --------------------------------------------------------------- the pipeline --

run_pipeline <- function(rec, cores, label = "pipeline") {
    bench(rec, label, "data(triticum)", utils::data("triticum",
        package = "MaMaMIA",
        envir = environment()
    ), size = FALSE)

    rca <- bench(
        rec, label, "RCA()",
        RCA(dcov = don_cov, dgc = don_gc, rcov = rec_cov, rgc = rec_gc, meta = meta)
    )

    corrected <- bench(
        rec, label, sprintf("correctReadCounts(cores=%d)", cores),
        correctReadCounts(rca, cores = cores, verbose = FALSE)
    )

    subset_cd <- bench(
        rec, label, "subset(At / A)",
        subset(corrected, don_subgenomes = "At", rec_subgenomes = "A")
    )

    bench(
        rec, label, "reverseWindows(RCA)",
        reverseWindows(subset_cd, chr_ids = "OY997261.1")
    )

    pairs <- stats::setNames(subset_cd$meta$rec_chr_ids, subset_cd$meta$don_chr_ids)
    isa <- bench(
        rec, label, "segments()",
        segments(subset_cd, target_pairs = pairs, verbose = FALSE)
    )

    bench(rec, label, "summary(ISA)", summary(isa, unit = "Mb", digits = 2))

    bedpe <- tempfile(fileext = ".bedpe")
    bench(rec, label, "writeToBedpe()", writeToBedpe(isa, file = bedpe), size = FALSE)

    bench(rec, label, "plot(ISA, mirror)",
        plot(isa, plot.type = "mirror"),
        size = FALSE
    )
    bench(rec, label, "plot(ISA, diff)",
        plot(isa, plot.type = "diff"),
        size = FALSE
    )

    list(rca = rca, corrected = corrected, subset_cd = subset_cd, isa = isa, bedpe = bedpe)
}

# Piecewise breakdown of correctReadCounts(), mirroring its body.
correction_internals <- function(rec, cores) {
    rca <- bench(rec, "correction", "data + RCA()", {
        utils::data("triticum", package = "MaMaMIA", envir = environment())
        RCA(dcov = don_cov, dgc = don_gc, rcov = rec_cov, rgc = rec_gc, meta = meta)
    }, size = FALSE)

    dat <- rca$data

    dat <- bench(rec, "correction", "  filter: valid/ideal + quantiles", {
        dat <- dplyr::mutate(dat, valid = (cov >= 0) & (gc > 0))
        gc_q <- stats::quantile(dat$gc[dat$valid], c(0.01, 0.99))
        cov_q <- stats::quantile(dat$cov[dat$valid], 0.99)
        dplyr::mutate(
            dat,
            ideal = valid & (gc >= gc_q[1]) & (gc <= gc_q[2]) & (cov <= cov_q)
        )
    })

    ideal_df <- dplyr::filter(dat, ideal)

    fit <- bench(
        rec, "correction", "  glmmTMB fit",
        glmmTMB::glmmTMB(
            cov ~ s(gc, k = 10) + subgenome,
            ziformula = ~ s(gc, k = 10) + subgenome,
            family = glmmTMB::nbinom2(),
            data = ideal_df, REML = TRUE,
            control = glmmTMB::glmmTMBControl(parallel = list(n = cores))
        )
    )

    gc_ref <- stats::median(dat$gc[dat$ideal], na.rm = TRUE)
    pred <- bench(
        rec, "correction", "  predict (fast path)",
        MaMaMIA:::predict_zinb_response(fit, dat$gc, dat$subgenome, gc_ref)
    )

    bench(rec, "correction", "  cor.gc + post-filter", {
        dat$cor.gc <- dat$cov * (pred$ref / (pred$actual + 1e-8))
        dat$ideal <- dat$ideal &
            dat$cor.gc < stats::quantile(dat$cor.gc, 0.99, na.rm = TRUE)
        dat
    })

    list(fit = fit, dat = dat, gc_ref = gc_ref, pred = pred)
}

# Detail of the prediction path, so its two costs are visible separately.
predict_detail <- function(rec, fit, dat, gc_ref) {
    key <- function(g, s) paste(sprintf("%.17g", g), s, sep = "\r")
    base <- data.frame(gc = dat$gc, subgenome = dat$subgenome)
    unique_rows <- base[!duplicated(key(base$gc, base$subgenome)), , drop = FALSE]
    ref_rows <- data.frame(gc = gc_ref, subgenome = unique(base$subgenome))
    newdata <- rbind(unique_rows, ref_rows)

    tmb <- bench(
        rec, "predict",
        sprintf("  predict(debug=TRUE) build (%d rows)", nrow(newdata)),
        stats::predict(fit, newdata = newdata, debug = TRUE),
        size = FALSE
    )

    bench(rec, "predict", "  design-matrix eta math", {
        pars <- fit$fit$parfull
        get_par <- function(name) pars[names(pars) == name]
        eta_cond <- as.numeric(as.matrix(tmb$data.tmb$X) %*% get_par("beta")) +
            as.numeric(as.matrix(tmb$data.tmb$Z) %*% get_par("b"))
        eta_zi <- as.numeric(as.matrix(tmb$data.tmb$Xzi) %*% get_par("betazi")) +
            as.numeric(as.matrix(tmb$data.tmb$Zzi) %*% get_par("bzi"))
        (1 - stats::plogis(eta_zi)) * exp(eta_cond)
    }, size = FALSE)

    bench(
        rec, "predict", "  glmmTMB::predict() (reference, both calls)",
        {
            nd <- data.frame(gc = dat$gc, subgenome = dat$subgenome)
            a <- stats::predict(fit, newdata = nd, type = "response")
            r <- stats::predict(fit, newdata = transform(nd, gc = gc_ref), type = "response")
            list(actual = as.numeric(a), ref = as.numeric(r))
        }
    )
}

# --------------------------------------------------------------- validation --

output_fingerprint <- function(res, cores) {
    segment_numbers <- res$isa$segments[vapply(
        res$isa$segments, is.numeric, logical(1L)
    )]
    c(
        cores = cores,
        n_ideal = sum(res$corrected$data$ideal),
        n_segments = nrow(res$isa$segments),
        n_putative = sum(res$isa$segments$putative_introgression),
        # Gating values: rounded to 10 significant digits. Rearranging the
        # arithmetic moves these outputs by up to 1.5e-14 relative (measured for
        # the direct-vs-glmmTMB design paths), so gating nearer the last bit only
        # catches values straddling a rounding boundary: that difference flipped
        # one of 35000 values at 12 digits and none at 11. Ten digits still
        # catches the nearest real signal, the ~1e-7 thread dependence.
        r_cor_gc = hash_object(signif(res$corrected$data$cor.gc, 10L)),
        r_fixef = hash_object(signif(unlist(glmmTMB::fixef(res$corrected$fit)), 10L)),
        r_segments = hash_object(lapply(segment_numbers, signif, 10L)),
        # Exact hashes, informational: these change whenever an optimisation
        # rearranges the arithmetic, even when no value changes appreciably.
        h_cor_gc = hash_object(res$corrected$data$cor.gc),
        h_ideal = hash_object(res$corrected$data$ideal),
        h_outliers = hash_object(res$corrected$outliers),
        h_fixef = hash_object(unlist(glmmTMB::fixef(res$corrected$fit))),
        h_segments = hash_object(res$isa$segments),
        h_bedpe = hash_object(readLines(res$bedpe))
    )
}

# Keys whose exact hash is reported but which are not allowed to fail --compare
# on their own: they are gated by their r_ counterpart instead.
INFORMATIONAL_KEYS <- c("h_cor_gc", "h_fixef", "h_segments")

read_section <- function(path, section) {
    lines <- readLines(path, warn = FALSE)
    start <- match(paste0("[", section, "]"), lines)
    if (is.na(start) || start == length(lines)) {
        return(character())
    }
    rest <- lines[(start + 1L):length(lines)]
    ends <- which(startsWith(rest, "["))
    if (length(ends)) {
        rest <- rest[seq_len(ends[1L] - 1L)]
    }
    rest <- rest[grepl(":", rest)]
    if (!length(rest)) {
        return(character())
    }
    stats::setNames(sub("^[^:]*:[[:space:]]*", "", rest), sub(":.*$", "", rest))
}

# Only the [outputs] section takes part in --compare: [env] legitimately differs
# between package versions and [checks] is derived from the run itself.
read_fingerprint <- function(path) read_section(path, "outputs")

# ---------------------------------------------------------------- reporting --

write_results <- function(rec, out_dir) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    csv <- file.path(out_dir, "RESULTS.csv")
    utils::write.csv(do.call(rbind, rec$rows), csv, row.names = FALSE)
    csv
}

write_fingerprint <- function(rec, fp, out_dir, opts) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    path <- file.path(out_dir, "FINGERPRINT.txt")
    env <- c(
        package_version = as.character(utils::packageVersion("MaMaMIA")),
        R_version = paste(R.version$major, R.version$minor, sep = "."),
        platform = R.version$platform,
        glmmTMB_version = as.character(utils::packageVersion("glmmTMB")),
        DNAcopy_version = as.character(utils::packageVersion("DNAcopy")),
        reps = opts$reps,
        seed = 1L
    )
    checks <- do.call(rbind, rec$checks)
    lines <- c(
        "# MaMaMIA pipeline fingerprint",
        "# Only the [outputs] section is compared by --compare; [env] is context.",
        "[env]",
        sprintf("%s: %s", names(env), env),
        "[outputs]",
        sprintf("%s: %s", names(fp), fp),
        "[checks]",
        sprintf("%s: %s", gsub(":", "_", checks$check), checks$result)
    )
    writeLines(lines, path)
    path
}

compare_fingerprints <- function(previous, current) {
    keys <- union(names(previous), names(current))
    do.call(rbind, lapply(keys, function(k) {
        p <- if (k %in% names(previous)) previous[[k]] else "<absent>"
        c <- if (k %in% names(current)) current[[k]] else "<absent>"
        data.frame(key = k, previous = p, current = c, match = identical(p, c))
    }))
}

# --------------------------------------------------------------------- main --

main <- function() {
    opts <- parse_args(commandArgs(trailingOnly = TRUE))
    suppressPackageStartupMessages(library(MaMaMIA))
    rec <- new_recorder()
    set.seed(1L)

    cat("== pipeline =====================================================\n")
    res <- run_pipeline(rec, opts$cores)
    print_section(rec, "pipeline")

    cat("\n== correctReadCounts() internals ================================\n")
    internals <- correction_internals(rec, opts$cores)
    print_section(rec, "correction")

    cat("\n== prediction path ==============================================\n")
    reference <- predict_detail(rec, internals$fit, internals$dat, internals$gc_ref)
    print_section(rec, "predict")

    cat("\n== consistency ==================================================\n")

    # 1. The fast path must agree with glmmTMB's own prediction.
    fast <- internals$pred
    max_rel <- max(
        abs(fast$actual - reference$actual) / pmax(abs(reference$actual), 1e-300)
    )
    check(
        rec, "fast prediction matches glmmTMB::predict()",
        max_rel < 1e-6,
        sprintf("max relative difference %.3e", max_rel)
    )

    # 2. The piecewise sum must account for the whole-call time.
    rows <- do.call(rbind, rec$rows)
    whole <- rows$elapsed_s[rows$step == sprintf("correctReadCounts(cores=%d)", opts$cores)]
    piecewise <- sum(rows$elapsed_s[
        rows$section == "correction" & grepl("^  ", rows$step)
    ])
    check(
        rec, "piecewise internals account for the whole call",
        abs(piecewise - whole) / whole < 0.15,
        sprintf("whole %.2f s vs piecewise %.2f s", whole, piecewise)
    )

    # 3. correctReadCounts() must be repeatable at fixed cores.
    repeat_fit <- suppressMessages(
        correctReadCounts(res$rca, cores = opts$cores, verbose = FALSE)
    )
    check(
        rec, "correctReadCounts() is repeatable",
        identical(
            hash_object(repeat_fit$data$cor.gc),
            hash_object(res$corrected$data$cor.gc)
        ),
        "identical corrected coverage for the same cores"
    )

    # 4. segments() must be reproducible and must not consume the caller's RNG.
    subset_cd <- res$subset_cd
    pairs <- stats::setNames(subset_cd$meta$rec_chr_ids, subset_cd$meta$don_chr_ids)
    set.seed(42L)
    first <- segments(subset_cd, target_pairs = pairs, verbose = FALSE)
    set.seed(99L)
    second <- segments(subset_cd, target_pairs = pairs, verbose = FALSE)
    check(
        rec, "segments() is independent of the caller's RNG state",
        identical(first$segments, second$segments),
        sprintf("%d segments from both calls", nrow(first$segments))
    )

    set.seed(123L)
    baseline <- stats::runif(3)
    set.seed(123L)
    invisible(segments(subset_cd, target_pairs = pairs, verbose = FALSE))
    after <- stats::runif(3)
    check(
        rec, "segments() leaves the RNG stream untouched",
        identical(baseline, after),
        "with the same seed, segments() does not advance the stream"
    )

    # 5. writeToBedpe() must be byte-stable.
    f1 <- tempfile(fileext = ".bedpe")
    f2 <- tempfile(fileext = ".bedpe")
    suppressMessages(writeToBedpe(res$isa, file = f1))
    suppressMessages(writeToBedpe(res$isa, file = f2))
    check(
        rec, "writeToBedpe() is byte-stable",
        identical(readLines(f1), readLines(f2)),
        sprintf("%d data rows", length(readLines(f1)) - 1L)
    )

    # 6. Whole-pipeline fingerprints must match across repeats.
    fingerprints <- list(output_fingerprint(res, opts$cores))
    if (opts$reps > 1L) {
        for (i in seq_len(opts$reps - 1L)) {
            set.seed(1L)
            cat(sprintf("   repeat %d of %d ...\n", i + 1L, opts$reps))
            again <- run_pipeline(rec, opts$cores, label = sprintf("repeat%d", i + 1L))
            fingerprints[[length(fingerprints) + 1L]] <-
                output_fingerprint(again, opts$cores)
        }
    }
    fp <- fingerprints[[1L]]
    same <- all(vapply(
        fingerprints[-1L],
        function(other) identical(other, fp),
        logical(1L)
    ))
    check(
        rec, sprintf("pipeline output identical across %d run(s)", opts$reps),
        same,
        "fingerprints compared"
    )

    if (!is.na(opts$compare)) {
        if (!file.exists(opts$compare)) abort("no such fingerprint: ", opts$compare)
        reference_env <- read_section(opts$compare, "env")
        cat(sprintf("\n== compared with %s ==\n", opts$compare))
        if (length(reference_env)) {
            cat(sprintf(
                "   reference env: %s\n",
                paste(sprintf("%s=%s", names(reference_env), reference_env), collapse = " ")
            ))
        }
        reference <- read_fingerprint(opts$compare)
        if (!length(reference)) {
            abort(
                "no [outputs] section in ", opts$compare,
                " -- was it written by an older version of this script?"
            )
        }
        missing <- setdiff(names(fp), names(reference))
        if (length(missing)) {
            abort(
                "reference fingerprint predates ", paste(missing, collapse = ", "),
                "; regenerate it with the current script"
            )
        }
        diff <- compare_fingerprints(reference, fp)
        gating <- !(diff$key %in% INFORMATIONAL_KEYS)
        n_diff <- sum(!diff$match)
        n_gate <- sum(!diff$match & gating)
        n_info <- sum(!diff$match & !gating)
        if (n_diff > 0L) {
            print(diff[!diff$match, , drop = FALSE], row.names = FALSE)
        }
        if (n_info > 0L && n_gate == 0L) {
            cat(sprintf(
                "   %d exact hash(es) differ but every value agrees to 10 significant digits\n",
                n_info
            ))
        }
        check(
            rec, "outputs unchanged versus reference fingerprint",
            n_gate == 0L,
            sprintf(
                "%d of %d gating keys differ, %d arithmetic-only",
                n_gate, sum(gating), n_info
            )
        )
    }

    print(do.call(rbind, rec$checks), row.names = FALSE)

    csv <- write_results(rec, opts$out)
    fp_path <- write_fingerprint(rec, fp, opts$out, opts)
    cat(sprintf("\nwrote %s\nwrote %s\n", csv, fp_path))

    failed <- vapply(rec$checks, function(x) x$result != "PASS", logical(1L))
    if (any(failed)) {
        cat("\nFAILED checks:\n")
        for (x in rec$checks[failed]) cat("  - ", x$check, ": ", x$detail, "\n", sep = "")
        quit(status = 1L)
    }
    cat("\nall consistency checks passed\n")
    invisible(NULL)
}

main()

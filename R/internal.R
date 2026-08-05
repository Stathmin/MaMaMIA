#' @noRd
validate_target_pairs <- function(target_pairs,
                                  don_chr_ids,
                                  rec_chr_ids) {
    if (length(target_pairs) < 1L || is.null(names(target_pairs)) ||
        any(!nzchar(names(target_pairs)))) {
        stop(
            "`target_pairs` must be a non-empty named vector of donor -> recipient chr_id",
            call. = FALSE
        )
    }
    
    target_don <- as.character(names(target_pairs))
    target_rec <- as.character(unname(target_pairs))
    stopifnot(
        "Some donor chromosome IDs in target_pairs are not present" =
            all(target_don %in% don_chr_ids)
    )
    stopifnot(
        "Some recipient chromosome IDs in target_pairs are not present" =
            all(target_rec %in% rec_chr_ids)
    )
    
    pairs <- data.frame(
        chr_id.don = target_don,
        chr_id.rec = target_rec,
        stringsAsFactors = FALSE
    )
    pairs <- pairs[!duplicated(pairs), , drop = FALSE]
    rownames(pairs) <- seq_len(nrow(pairs))
    pairs
}

#' @noRd
empty_segments <- function() {
    data.frame(
        chr_id.don = character(),
        chr_id.rec = character(),
        loc.start = integer(),
        loc.end = integer(),
        num.mark = integer(),
        seg.mean = numeric(),
        seg.sd = numeric(),
        seg.median = numeric(),
        seg.mad = numeric(),
        putative_introgression = logical(),
        chr_name.don = character(),
        subgenome.don = character(),
        chr_name.rec = character(),
        subgenome.rec = character(),
        stringsAsFactors = FALSE
    )
}

#' @noRd
annotate_segment_coords <- function(segments, out) {
    merge_data <- dplyr::select(out, dplyr::all_of(
        c(
            'iid',
            'chr_id.don',
            'start.don',
            'stop.don',
            'chr_id.rec',
            'start.rec',
            'stop.rec'
        )
    ))
    segments |>
        dplyr::inner_join(
            y = dplyr::select(merge_data, !dplyr::all_of(c(
                'stop.don', 'stop.rec'
            ))),
            by = c('chr_id.don', 'chr_id.rec', 'loc.start' = 'iid')
        ) |>
        dplyr::inner_join(
            y = dplyr::select(merge_data, !dplyr::all_of(c(
                'start.don', 'start.rec'
            ))),
            by = c('chr_id.don', 'chr_id.rec', 'loc.end' = 'iid')
        )
}

#' @noRd
getSegments <- function(counts, chrom, maploc, alpha, undo.SD) {
    CNA.object <- DNAcopy::CNA(
        genomdat = counts,
        chrom = chrom,
        maploc = maploc,
        data.type = 'logratio'
    ) |> DNAcopy::smooth.CNA()
    DNAcopy::segment(
        CNA.object,
        verbose = 0,
        undo.splits = 'sdundo',
        alpha = alpha,
        min.width = 2,
        undo.SD = undo.SD
    )
}

#' @noRd
validate_segmentation_params <- function(alpha, min_width, undo_SD) {
    if (alpha <= 0 || alpha > 1) {
        stop('1 >= alpha > 0 is not satisfied', call. = FALSE)
    }
    
    if (min_width < 2) {
        stop('min_width must be >=2', call. = FALSE)
    }
    
    if (undo_SD <= 0 || undo_SD > 10) {
        stop('10 >= undo_SD > 0 is not satisfied', call. = FALSE)
    }
    
    list(alpha = alpha, undo.SD = undo_SD)
}


#' @noRd
as_window_df <- function(x, value_col, arg_name) {
    if (inherits(x, "matrix")) {
        x <- as.data.frame(x, stringsAsFactors = FALSE)
    }
    if (!is.data.frame(x)) {
        stop(sprintf("`%s` must be a data.frame, tibble, or matrix", arg_name),
             call. = FALSE)
    }
    if (ncol(x) != 4L) {
        stop(
            sprintf(
                "`%s` must have exactly 4 columns (chr_id, start, stop, %s)",
                arg_name,
                value_col
            ),
            call. = FALSE
        )
    }
    
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    names(x) <- c("chr_id", "start", "stop", value_col)
    
    if (is.factor(x$chr_id))
        x$chr_id <- as.character(x$chr_id)
    if (!is.atomic(x$chr_id) || is.list(x$chr_id)) {
        stop(
            sprintf(
                "`%s$chr_id` must be atomic (character/integer/factor)",
                arg_name
            ),
            call. = FALSE
        )
    }
    x$chr_id <- as.character(x$chr_id)
    
    for (col in c("start", "stop", value_col)) {
        if (!is.numeric(x[[col]])) {
            coerced <- suppressWarnings(as.numeric(x[[col]]))
            if (anyNA(coerced) && !anyNA(x[[col]])) {
                stop(sprintf("`%s$%s` must be numeric", arg_name, col),
                     call. = FALSE)
            }
            x[[col]] <- coerced
        }
    }
    
    if (anyNA(x$start) || anyNA(x$stop)) {
        stop(sprintf("`%s` start/stop must not contain NA", arg_name),
             call. = FALSE)
    }
    if (anyNA(x[[value_col]])) {
        stop(
            sprintf("`%s$%s` must not contain NA", arg_name, value_col),
            call. = FALSE
        )
    }
    bad <- which(!(x$start < x$stop))
    if (length(bad) > 0L) {
        stop(
            sprintf(
                "`%s` has intervals with start >= stop (first at row %d)",
                arg_name,
                bad[[1L]]
            ),
            call. = FALSE
        )
    }
    
    return(x)
}

#' @noRd
validate_windows_tiled <- function(x, arg_name) {
    for (chr in unique(x$chr_id)) {
        idx <- which(x$chr_id == chr)
        if (length(idx) < 2L) {
            next
        }
        start <- x$start[idx]
        stop <- x$stop[idx]
        next_start <- start[-1L]
        prev_stop <- stop[-length(stop)]
        if (any(next_start < prev_stop)) {
            stop(
                sprintf(
                    "`%s` has overlapping windows on chromosome %s",
                    arg_name,
                    chr
                ),
                call. = FALSE
            )
        }
        if (any(next_start > prev_stop)) {
            stop(
                sprintf(
                    paste0(
                        "`%s` has gaps between windows on chromosome %s ",
                        "(expected contiguous half-open intervals)"
                    ),
                    arg_name,
                    chr
                ),
                call. = FALSE
            )
        }
    }
    invisible(x)
}


#' @noRd
as_meta_df <- function(x, arg_name = "meta") {
    if (inherits(x, "matrix")) {
        x <- as.data.frame(x, stringsAsFactors = FALSE)
    }
    if (!is.data.frame(x)) {
        stop(sprintf("`%s` must be a data.frame, tibble, or matrix", arg_name),
             call. = FALSE)
    }
    if (ncol(x) != 3L) {
        stop(
            sprintf(
                "`%s` must have exactly 3 columns (chr_id, chr_name, subgenome)",
                arg_name
            ),
            call. = FALSE
        )
    }
    
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    names(x) <- c("chr_id", "chr_name", "subgenome")
    
    for (col in c("chr_id", "chr_name", "subgenome")) {
        if (is.factor(x[[col]]))
            x[[col]] <- as.character(x[[col]])
        if (!is.atomic(x[[col]]) || is.list(x[[col]])) {
            stop(
                sprintf(
                    "`%s$%s` must be atomic (character/integer/factor)",
                    arg_name,
                    col
                ),
                call. = FALSE
            )
        }
        x[[col]] <- as.character(x[[col]])
    }
    
    if (anyNA(x)) {
        stop(sprintf("`%s` must not contain NA", arg_name), call. = FALSE)
    }
    if (anyDuplicated(x$chr_id)) {
        stop(sprintf("`%s$chr_id` must be unique", arg_name), call. = FALSE)
    }
    
    return(x)
}

#' @noRd
resolve_chr_targets <- function(meta,
                                don_chromlist,
                                rec_chromlist,
                                don_subgenomes,
                                rec_subgenomes,
                                don_chr_ids,
                                rec_chr_ids) {
    stopifnot(
        "Require at least one of don_chromlist, rec_chromlist, don_subgenomes, rec_subgenomes" = !all(vapply(
            list(
                don_chromlist,
                rec_chromlist,
                don_subgenomes,
                rec_subgenomes
            ),
            is.null,
            logical(1)
        ))
    )
    
    if (!is.null(don_chromlist) && length(don_chromlist) == 0L) {
        stop("`don_chromlist` must be NULL or non-empty", call. = FALSE)
    }
    if (!is.null(rec_chromlist) && length(rec_chromlist) == 0L) {
        stop("`rec_chromlist` must be NULL or non-empty", call. = FALSE)
    }
    if (!is.null(don_subgenomes) && length(don_subgenomes) == 0L) {
        stop("`don_subgenomes` must be NULL or non-empty", call. = FALSE)
    }
    if (!is.null(rec_subgenomes) && length(rec_subgenomes) == 0L) {
        stop("`rec_subgenomes` must be NULL or non-empty", call. = FALSE)
    }
    
    if (!is.null(don_subgenomes)) {
        stopifnot("Unknown donor subgenome label(s)" =
                      all(don_subgenomes %in% meta$subgenome))
    }
    if (!is.null(rec_subgenomes)) {
        stopifnot("Unknown recipient subgenome label(s)" =
                      all(rec_subgenomes %in% meta$subgenome))
    }
    if (!is.null(don_chromlist)) {
        stopifnot("Unknown donor chromosome ID(s)" =
                      all(don_chromlist %in% don_chr_ids))
    }
    if (!is.null(rec_chromlist)) {
        stopifnot("Unknown recipient chromosome ID(s)" =
                      all(rec_chromlist %in% rec_chr_ids))
    }
    
    resolve_one <- function(chromlist, subgenomes, current_ids) {
        current_ids <- as.character(current_ids)
        if (is.null(chromlist) && is.null(subgenomes)) {
            return(current_ids)
        }
        ids <- character(0)
        if (!is.null(subgenomes)) {
            ids <- c(
                ids,
                intersect(
                    as.character(meta$chr_id[meta$subgenome %in% subgenomes]),
                    current_ids
                )
            )
        }
        if (!is.null(chromlist)) {
            ids <- c(ids, as.character(chromlist))
        }
        unique(ids)
    }
    
    target_don <- resolve_one(don_chromlist, don_subgenomes, don_chr_ids)
    target_rec <- resolve_one(rec_chromlist, rec_subgenomes, rec_chr_ids)
    if (length(target_don) == 0L || length(target_rec) == 0L) {
        stop(
            "Subset would leave donor or recipient with no chromosomes",
            call. = FALSE
        )
    }
    
    list(target_don, target_rec)
}


#' @noRd
validate_equal_n_windows <- function(chr_id) {
    counts <- table(chr_id, useNA = "no")
    if (length(counts) == 0L) {
        stop("Input must contain at least one chromosome", call. = FALSE)
    }
    n_unique <- length(unique(as.integer(counts)))
    if (n_unique != 1L) {
        detail <- paste(
            sprintf("%s=%d", names(counts), as.integer(counts)),
            collapse = ", "
        )
        stop(
            paste0(
                "All chromosomes must have the same number of windows; found: ",
                detail
            ),
            call. = FALSE
        )
    }
    invisible(as.integer(counts[[1L]]))
}

#' @noRd
check_cores <- function(x) {
    x <- as.integer(x)[1L]
    max_cores <- parallel::detectCores()
    if (is.na(max_cores) || max_cores < 1L) {
        max_cores <- 1L
    }
    if (is.na(x) || x < 1L || x > max_cores) {
        warning(
            paste0(
                'Cores must be between 1 and ',
                max_cores,
                '. Fallback to 1 core.'
            ),
            call. = FALSE
        )
        return(1L)
    }
    return(x)
}


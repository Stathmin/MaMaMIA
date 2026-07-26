get_segments <- function(counts, chrom, maploc, param) {
    CNA.object <- DNAcopy::CNA(
        genomdat = counts,
        chrom = chrom,
        maploc = maploc,
        data.type = 'logratio'
    )
    DNAcopy::segment(
        CNA.object,
        verbose = 0,
        undo.splits = 'sdundo',
        alpha = param$alpha,
        min.width = pmin(5, param$min.width),
        undo.SD = param$undo.SD
    )
}


segments <- function(RCA,
                     param = NULL,
                     getparam = FALSE,
                     verbose = TRUE) {
    if (getparam == FALSE) {
        param <- list(alpha = 0.001, # Less value - less segments, but may miss small introgressions
                      min.width = 5, # The minimum number of windows in a changed segment
                      undo.SD = 1 # The number of standard deviations between segment means required to keep a split)
                      )}
    
    if (verbose == TRUE) {
        message('Performing segmentation...')
    }
    prepared_data <- tidyr::expand_grid(
        chr_id_don = RCA$meta$don_chr_ids,
        chr_id_rec = RCA$meta$rec_chr_ids
    ) |>
        dplyr::mutate(
            data = purrr::map2(chr_id_don, chr_id_rec, \(x, y) {
                don_data <- dplyr::filter(RCA$data, chr_id == x) |>
                    dplyr::rename_with( ~ paste0(., ".don"))
                rec_data <- dplyr::filter(RCA$data, chr_id == y) |>
                    dplyr::rename_with( ~ paste0(., ".rec"))
                dplyr::bind_cols(don_data, rec_data) |>
                    dplyr::mutate(diff = cor.gc.don - cor.gc.rec) |>
                    dplyr::filter(ideal.don & ideal.rec)
            }),
            segments = purrr::map(data, \(x) {
                get_segments(
                    counts = x$diff,
                    chrom = paste(x$chr_id.don, x$chr_id.rec, sep = '>'),
                    maploc = x$iid.don,
                    param = param
                ) |>
                    DNAcopy::segments.summary() |>
                    dplyr::select(-ID) |>
                    dplyr::filter(seg.median >= 0) |>
                    tidyr::separate(chrom, c('chr_id.don', 'chr_id.rec'), sep = '>') |>
                    dplyr::inner_join(RCA$meta$meta, by = c('chr_id.don' = 'chr_id')) |>
                    dplyr::inner_join(
                        RCA$meta$meta,
                        by = c('chr_id.rec' = 'chr_id'),
                        suffix = c('.don', '.rec')
                    )
            }),
            .keep = 'unused'
        )
    
    structure(
        list(
            data = RCA$data,
            out = dplyr::bind_rows(prepared_data$data),
            segments = dplyr::bind_rows(prepared_data$segments),
            meta = RCA$meta
        ),
        class = 'ISA'
    )
}


subset.ISA <- function(ISA,
                       don_chromlist = NULL,
                       rec_chromlist = NULL,
                       don_subgenomes = NULL,
                       rec_subgenomes = NULL) {
    if (!is.null(don_subgenomes) | !is.null(rec_subgenomes)) {
        target_don_chrs <- c(ISA$meta$meta$chr_id[ISA$meta$meta$subgenome %in% don_subgenomes], don_chromlist) |>
            unique()
        target_rec_chrs <- c(ISA$meta$meta$chr_id[ISA$meta$meta$subgenome %in% rec_subgenomes], rec_chromlist) |>
            unique()
    } else {
        target_don_chrs <- unique(don_chromlist)
        target_rec_chrs <- unique(rec_chromlist)
    }
    
    ISA$data <- dplyr::filter(ISA$data, chr_id %in% c(target_don_chrs, target_rec_chrs))
    ISA$out <- dplyr::filter(ISA$out,
                             chr_id.don %in% target_don_chrs &
                                 chr_id.rec %in% target_rec_chrs)
    ISA$segments <- dplyr::filter(ISA$segments,
                                  chr_id.don %in% target_don_chrs &
                                      chr_id.rec %in% target_rec_chrs)
    ISA$meta$meta <- dplyr::filter(ISA$meta$meta,
                                   chr_id %in% c(target_don_chrs, target_rec_chrs))
    ISA$meta$rec_chr_ids <- target_rec_chrs
    ISA$meta$don_chr_ids <- target_don_chrs
    
    return(ISA)
}


plot.ISA <- function(ISA,
                     plot.type = c('mirror', 'diff'),
                     target_pairs = NULL,
                     ...) {
    plot_type <- match.arg(plot.type)
    
    if (!is.null(target_pairs)) {
        ISA$out <- dplyr::semi_join(
            ISA$out,
            data.frame(chr_id.don = names(target_pairs), chr_id.rec = target_pairs),
            by = c('chr_id.don', 'chr_id.rec')
        )
        
        ISA$segments <- dplyr::semi_join(
            ISA$segments,
            data.frame(chr_id.don = names(target_pairs), chr_id.rec = target_pairs),
            by = c('chr_id.don', 'chr_id.rec')
        )
    }
    
    # Add throughout number of windows
    ISA$out$iid.e2e <- seq_along(ISA$out$iid.don)
    
    # Convert windows coordinates into end-to-end scale
    ISA$segments <- dplyr::inner_join(
        ISA$segments,
        dplyr::select(ISA$out, chr_id.don, chr_id.rec, iid.don, iid.e2e),
        by = c(
            'chr_id.don' = 'chr_id.don',
            'chr_id.rec' = 'chr_id.rec',
            'loc.start' = 'iid.don'
        )
    ) |>
        dplyr::inner_join(
            dplyr::select(ISA$out, chr_id.don, chr_id.rec, iid.don, iid.e2e),
            by = c(
                'chr_id.don' = 'chr_id.don',
                'chr_id.rec' = 'chr_id.rec',
                'loc.end' = 'iid.don'
            ),
            suffix = c('.start', '.end')
        )
    
    # Color of each chromosome
    ISA$out$color.don <- (match(ISA$out$chr_name.don, unique(ISA$out$chr_name.don)) - 1) %% 2
    ISA$out$color.rec <- (match(ISA$out$chr_name.don, unique(ISA$out$chr_name.don)) - 1) %% 2
    
    # Calculate break coordinates
    breaks <- ISA$out |>
        dplyr::group_by(chr_name.don, chr_name.rec) |>
        dplyr::reframe(pos = min(iid.e2e) + (diff(range(iid.e2e)) / 2))
    
    plot <- ggplot2::ggplot(data = ISA$out) + 
        ggplot2::theme_bw() +
        ggplot2::scale_x_continuous(
            expand = ggplot2::expansion(mult = c(0.01, 0.01)),
            breaks = breaks$pos,
            labels = paste(breaks$chr_name.don, breaks$chr_name.rec, sep = '/')
        ) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = c(1, 1)),
                                    labels = \(x) gsub('-', '', x),
                                    name = 'Coverage') +
        ggplot2::theme(aspect.ratio = 1 / 5,
                       axis.title.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank(),
                       panel.grid = ggplot2::element_blank(),
                       legend.position = 'none')

    if (plot_type == 'mirror') {
        plot + ggplot2::geom_point(mapping = ggplot2::aes(
            x = iid.e2e,
            y = cor.gc.don + 1,
            color = factor(color.don)
        ),
        size = 0.5) +
            ggplot2::scale_color_manual(values = c("grey60", "grey80")) +
            ggnewscale::new_scale_color() +
            ggplot2::geom_point(mapping = ggplot2::aes(
                x = iid.e2e,
                y = -cor.gc.rec - 1,
                color = factor(color.rec)
            ),
            size = 0.5) + 
            ggplot2::scale_color_manual(values = c("grey80", "grey60"))
    } else {
        plot + ggplot2::geom_point(mapping = ggplot2::aes(
            x = iid.e2e,
            y = diff,
            color = factor(color.don)
        ),
        size = 0.5) +
            ggplot2::geom_hline(yintercept = 0,
                                color = 'green',
                                linetype = 'dashed',
                                linewidth = 1) +
            ggplot2::geom_segment(
                data = ISA$segments,
                mapping = ggplot2::aes(
                    x = iid.e2e.start,
                    y = seg.median,
                    xend = iid.e2e.end,
                    yend = seg.median
                ),
                color = 'orangered',
                linewidth = 1
            ) +
            ggplot2::scale_color_manual(values = c("grey60", "grey80"))
    }
}

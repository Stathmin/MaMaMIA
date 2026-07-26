RCA <- function(dcov, dgc, rcov, rgc, meta, presorted = TRUE) {
    # RCA - Read Count Array
    
    # Same dims
    stopifnot("dcov and dgc has different number of columns and row" = dim(dcov) == dim(dgc))
    stopifnot("rcov and rgc has different number of columns and row" = dim(rcov) == dim(rgc))
    
    # Check if 4 columns are presented
    stopifnot("Coverage data must contain exactly 4 columns" = ncol(dcov) == 4)
    stopifnot("GC data must contain exactly 4 columns" = ncol(dgc) == 4)
    
    dcov <- `colnames<-`(dcov, c('chr_id', 'start', 'stop', 'cov'))
    rcov <- `colnames<-`(rcov, c('chr_id', 'start', 'stop', 'cov'))
    dgc <- `colnames<-`(dgc, c('chr_id', 'start', 'stop', 'gc'))
    rgc <- `colnames<-`(rgc, c('chr_id', 'start', 'stop', 'gc'))
    
    if (!presorted) {
        message('Sotring...')
        dcov <- dplyr::arrange(dcov, chr_id, start, stop)
        rcov <- dplyr::arrange(rcov, chr_id, start, stop)
        dgc <- dplyr::arrange(dgc, chr_id, start, stop)
        rgc <- dplyr::arrange(rgc, chr_id, start, stop)
    }
    
    # Same chromosomes
    stopifnot("Chromosome IDs are not matching in donor coverage and GC data" = all(dcov[, 1] == dgc[, 1]))
    stopifnot("Chromosome IDs are not matching in recipient coverage and GC data" = all(rcov[, 1] == rgc[, 1]))
    
    # Same coordinates
    stopifnot("Window coordinates are not matching in donor coverage and GC data" = all(dcov[, 2:3] == dgc[, 2:3]))
    stopifnot("Window coordinates are not matching in donor coverage and GC data" = all(rcov[, 2:3] == rgc[, 2:3]))
    
    # Validate meta info
    stopifnot("Metainfo must not contain NA" = sum(is.na(meta)) == 0)
    stopifnot("Metainfo must contain exactly 3 columns" = ncol(meta) == 3)
    stopifnot(
        "Column names of meta do not match specification" =
            colnames(meta) == c("chr_id", "chr_name", "subgenome")
    )
    stopifnot("Duplicated chromosomes in meta" = length(meta[['chr_id']]) == length(unique(meta[['chr_id']])))
    
    don_chr_ids <- unique(dcov$chr_id)
    rec_chr_ids <- unique(rcov$chr_id)
    stopifnot(
        "Number of chromosomes in meta do not match sum of donor and recipient chromosomes" =
            length(meta[['chr_id']]) == length(c(don_chr_ids, rec_chr_ids))
    )
    
    
    don_data <- dplyr::inner_join(dcov, dgc, by = c('chr_id', 'start', 'stop'))
    rec_data <- dplyr::inner_join(rcov, rgc, by = c('chr_id', 'start', 'stop'))
    full_data <- dplyr::bind_rows(rec_data, don_data) |>
        dplyr::full_join(meta, by = 'chr_id') |> 
        dplyr::group_by(chr_id) |> 
        dplyr::mutate(iid = dplyr::row_number(), .after = 'chr_id') |> 
        dplyr::ungroup()
    
    structure(list(
        data = full_data,
        meta = list(
            meta = meta,
            don_chr_ids = don_chr_ids,
            rec_chr_ids = rec_chr_ids
        )
    ), class = 'RCA')
}


correct_read_counts <- function(RCA,
                                cov_outlier = 0.01,
                                gc_outlier = 0.01,
                                verbose = TRUE) {
    if (verbose) {
        message('Applying filter on data...')
    }
    RCA$data <- RCA$data |>
        dplyr::mutate(valid = (cov >= 0) & (gc > 0)) |>
        dplyr::mutate(ideal = valid &
                          (gc >= quantile(gc[valid], gc_outlier)) &
                          (gc <= quantile(gc[valid], 1 - gc_outlier)) &
                          (cov <= quantile(cov[valid], 1 - cov_outlier)))
    
    
    if (verbose) {
        message('Correcting for GC bias...')
    }
    fit <- glmmTMB::glmmTMB(
        cov ~ s(gc + (1 | subgenome), k = 10),
        ziformula = ~ s(gc + (1 | subgenome), k = 10),
        family = glmmTMB::nbinom2(),
        data = dplyr::filter(RCA$data, ideal),
        REML = TRUE
    )
    
    mean_gc_ideal <- dplyr::filter(RCA$data, ideal) |>
        dplyr::pull(gc) |>
        mean()
    
    cor.gc <- predict(
        object = fit,
        newdata = transform(RCA$data, covariate = mean_gc_ideal),
        type = 'response'
    )
    
    RCA$data$cor.gc <- RCA$data$cov / cor.gc
    RCA$data$ideal <- RCA$data$ideal &
        RCA$data$cor.gc < quantile(RCA$data$cor.gc,
                                           probs = 1 - cov_outlier,
                                           na.rm = TRUE)
    return(RCA)
}


plot.RCA <- function(RCA,
                     plot.type = c('orig_cov', 'corr_cov'),
                     show_outliers = FALSE,
                     ...) {
    if (show_outliers == TRUE) {
        plotting_df <- RCA$data
    } else {
        plotting_df <- dplyr::filter(RCA$data, ideal == TRUE)
    }
    plotting_df$color <- (match(plotting_df$chr_id, unique(plotting_df$chr_id)) - 1) %% 2
    
    plot_type <- match.arg(plot.type)
    
    breaks <- plotting_df |>
        dplyr::mutate(rn = seq_along(chr_name)) |>
        dplyr::group_by(chr_name) |>
        dplyr::reframe(pos = min(rn) + (diff(range(rn)) / 2))
    
    plot <- ggplot2::ggplot(plotting_df, ggplot2::aes(x = seq_along(chr_id), color = factor(color))) +
        ggplot2::theme_bw() +
        ggplot2::scale_x_continuous(
            expand = ggplot2::expansion(),
            breaks = breaks$pos,
            labels = breaks$chr_name
        ) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = c(0, 1))) +
        ggplot2::scale_color_manual(values = c("grey60", "grey80")) +
        ggplot2::theme(
            aspect.ratio = 1 / 5,
            legend.position = 'none',
            panel.grid = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank()
        )
    
    if (plot_type == 'orig_cov') {
        plot + ggplot2::geom_point(size = 0.5, mapping = ggplot2::aes(y = cov))
    } else {
        plot + ggplot2::geom_point(size = 0.5, mapping = ggplot2::aes(y = cor.gc))
    }
}


subset.RCA <- function(RCA,
                       don_chromlist = NULL,
                       rec_chromlist = NULL,
                       don_subgenomes = NULL,
                       rec_subgenomes = NULL) {
    if (!is.null(don_subgenomes) | !is.null(rec_subgenomes)) {
        target_don_chrs <- c(RCA$meta$meta$chr_id[RCA$meta$meta$subgenome %in% don_subgenomes], don_chromlist) |>
            unique()
        target_rec_chrs <- c(RCA$meta$meta$chr_id[RCA$meta$meta$subgenome %in% rec_subgenomes], rec_chromlist) |>
            unique()
    } else {
        target_don_chrs <- unique(don_chromlist)
        target_rec_chrs <- unique(rec_chromlist)
    }
    
    RCA$data <- dplyr::filter(RCA$data, chr_id %in% c(target_don_chrs, target_rec_chrs))
    RCA$meta$meta <- dplyr::filter(RCA$meta$meta,
                                   chr_id %in% c(target_don_chrs, target_rec_chrs))
    RCA$meta$rec_chr_ids <- target_rec_chrs
    RCA$meta$don_chr_ids <- target_don_chrs
    
    return(RCA)
}

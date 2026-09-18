# Render optional vignette figures into vignettes/figures/.
# Run from package root; use --install for plotting deps, --rebuild-data to refit.

cran_packages <- c(
    "ggplot2", "dplyr", "ggtext", "gggenes", "ggdist", "ggraph",
    "tidygraph", "circlize", "patchwork", "png", "ragg"
)
args <- commandArgs(trailingOnly = TRUE)
if ("--install" %in% args) {
    missing <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing) > 0L) {
        install.packages(missing, repos = "https://packagemanager.posit.co/cran/latest")
    }
}
missing <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
    stop(
        "Missing illustration packages: ", paste(missing, collapse = ", "),
        "\nInstall them with: Rscript tools/illustrations/render-illustrations.R --install",
        call. = FALSE
    )
}

suppressMessages({
    library(ggplot2)
    library(dplyr)
    library(ggtext)
    library(gggenes)
    library(ggdist)
    library(ggraph)
    library(tidygraph)
    library(circlize)
    library(patchwork)
    library(png)
    library(ragg)
})
options(bitmapType = "cairo")

OUT <- "vignettes/figures"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
DATA_CACHE <- "/tmp/mamamia-illus.rds"
FIT_CACHE <- "/tmp/mamamia-fit.rds"
PAL <- list(don = "#d6604d", rec = "#4393c3", ink = "#333333",
            soft = "#9e9e9e", bg = "#f0f0f0")

theme_illus <- function(base = 11) {
    theme_minimal(base_size = base) +
        theme(
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(colour = "grey93"),
            plot.title = ggtext::element_markdown(face = "bold", size = base + 2),
            plot.subtitle = ggtext::element_markdown(colour = "grey30"),
            plot.caption = ggtext::element_markdown(
                colour = "grey40", hjust = 0, size = base - 2
            ),
            legend.position = "bottom"
        )
}

save_fig <- function(plot, name, width, height) {
    f <- file.path(OUT, name)
    ggsave(f, plot,
        width = width, height = height, dpi = 200, bg = "white",
        device = ragg::agg_png
    )
    cat("wrote", f, "\n")
}

# ---------------------------------------------------------------- data ------

build_illustration_data <- function(path = DATA_CACHE, fit_cache = FIT_CACHE) {
    suppressMessages(library(MaMaMIA))
    data("triticum", package = "MaMaMIA", envir = environment())
    rca <- RCA(
        dcov = don_cov, dgc = don_gc,
        rcov = rec_cov, rgc = rec_gc,
        meta = meta
    )
    cd <- correctReadCounts(
        rca, cores = 4L, verbose = FALSE,
        cache = fit_cache
    )
    sc <- subset(cd, don_subgenomes = "At", rec_subgenomes = "A")
    tp <- stats::setNames(sc$meta$rec_chr_ids, sc$meta$don_chr_ids)
    isa <- segments(sc, target_pairs = tp, verbose = FALSE)
    saveRDS(list(cd = cd, sc = sc, isa = isa), path)
    cat(
        "saved ", path, ": ", sum(cd$data$ideal), " ideal windows, ",
        nrow(isa$segments), " segments, ",
        sum(isa$segments$putative_introgression), " putative introgressions\n",
        sep = ""
    )
}

if (!file.exists(DATA_CACHE) || "--rebuild-data" %in% args) {
    build_illustration_data()
}
x <- readRDS(DATA_CACHE)
cd <- x$cd
sc <- x$sc
isa <- x$isa

chrs <- unique(sc$data[, c("chr_id", "chr_name")])
chr_len <- tapply(sc$data$stop, sc$data$chr_id, max)
chrs$len <- as.numeric(chr_len[as.character(chrs$chr_id)])
pairs <- isa$meta$target_pairs
don_names <- chrs$chr_name[match(pairs$chr_id.don, chrs$chr_id)]
rec_names <- chrs$chr_name[match(pairs$chr_id.rec, chrs$chr_id)]

pair_out <- function(i) {
    isa$out[isa$out$chr_id.don == pairs$chr_id.don[i] &
        isa$out$chr_id.rec == pairs$chr_id.rec[i], , drop = FALSE]
}
pair_seg <- function(i) {
    isa$segments[isa$segments$chr_id.don == pairs$chr_id.don[i] &
        isa$segments$chr_id.rec == pairs$chr_id.rec[i], , drop = FALSE]
}

# ------------------------------------------------- 1. how the signal arises --

fig_mechanism <- function() {
    out <- pair_out(1)
    seg <- pair_seg(1)
    called <- seg[seg$putative_introgression, , drop = FALSE]
    x_min <- min(out$iid)
    x_max <- max(out$iid)
    lo <- min(called$loc.start)
    hi <- max(called$loc.end)
    caveat <- numeric(0)

    schema_df <- data.frame(
        xmin = c(x_min, lo),
        xmax = c(x_max, hi),
        y = 1,
        forward = TRUE,
        type = c("recipient background", "introgressed segment"),
        stringsAsFactors = FALSE
    )
    schema <- ggplot(
        schema_df,
        aes(
            xmin = .data$xmin, xmax = .data$xmax, y = .data$y,
            forward = .data$forward, fill = .data$type
        )
    ) +
        geom_gene_arrow(
            arrowhead_height = unit(7, "mm"), arrowhead_width = unit(7, "mm")
        ) +
        scale_fill_manual(
            name = NULL,
            values = c(
                "recipient background" = "grey85",
                "introgressed segment" = PAL$don
            )
        ) +
        annotate("text",
            x = x_min + (x_max - x_min) * 0.30, y = 1.00, label = "recipient background",
            size = 3, colour = "grey25", vjust = 0.5
        ) +
        annotate("text",
            x = (lo + hi) / 2, y = 1.62,
            label = "introgressed segment", size = 3.2, colour = PAL$don
        ) +
        annotate("segment",
            x = (lo + hi) / 2, xend = (lo + hi) / 2,
            y = 1.50, yend = 1.22, colour = PAL$don,
            arrow = arrow(length = unit(0.10, "cm"))
        ) +
        scale_x_continuous(expand = c(0, 0), limits = c(x_min, x_max)) +
        scale_y_continuous(limits = c(0.6, 1.9)) +
        theme_genes() +
        theme(
            legend.position = "none",
            axis.text = element_blank(),
            axis.title = element_blank(),
            plot.margin = margin(0, 12, 0, 12)
        )

    cov <- data.frame(
        x = out$iid, donor = out$cor.gc.don, recipient = -out$cor.gc.rec
    )
    ymax <- max(abs(c(cov$donor, cov$recipient))) * 1.08

    profile <- ggplot(cov) +
        geom_rect(
            data = data.frame(xmin = lo, xmax = hi),
            aes(xmin = .data$xmin, xmax = .data$xmax),
            ymin = -ymax, ymax = ymax, fill = PAL$don, alpha = 0.10,
            inherit.aes = FALSE
        ) +
        geom_area(aes(x = .data$x, y = .data$donor),
            fill = PAL$don, alpha = 0.55
        ) +
        geom_area(aes(x = .data$x, y = .data$recipient),
            fill = PAL$rec, alpha = 0.55
        ) +
        geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.3) +
        annotate("text",
            x = x_min + (x_max - x_min) * 0.02,
            y = ymax * 0.78, label = "donor coverage",
            hjust = 0, size = 3.2, colour = PAL$don
        ) +
        annotate("text",
            x = x_min + (x_max - x_min) * 0.02,
            y = -ymax * 0.78, label = "recipient coverage",
            hjust = 0, size = 3.2, colour = PAL$rec
        ) +
        annotate("text",
            x = (lo + hi) / 2, y = ymax * 0.88,
            label = "donor above recipient:\nthis interval is called",
            size = 2.9, colour = PAL$ink, lineheight = 1.1
        ) +
        scale_x_continuous(expand = c(0, 0), limits = c(x_min, x_max)) +
        scale_y_continuous(
            limits = c(-ymax, ymax),
            labels = function(v) round(abs(v))
        ) +
        labs(
            x = "Window along the chromosome",
            y = "GC-corrected coverage\n(donor above, recipient below)"
        ) +
        theme_illus()

    (schema / profile) +
        plot_layout(heights = c(1, 4.2)) +
        plot_annotation(
            title = paste0(
                "An introgression shows up as donor excess: *", don_names[1],
                "* (donor) against *", rec_names[1], "* (recipient)"
            )
        ) &
        theme(
            plot.title = ggtext::element_markdown(face = "bold", size = 13),
            plot.margin = margin(6, 12, 6, 8)
        )
}

# ------------------------------------------------------------------ 2. flow --

fig_workflow <- function() {
    nodes <- data.frame(
        name = seq_len(7),
        label = c(
            "Sample\nrecipient chromosome carrying\na donor-derived segment",
            "Reads mapped competitively\nto a concatenated donor + recipient reference",
            "Windowed tables\nread counts and GC content, one table per genome",
            "GC-bias correction\nzero-inflated negative binomial model, giving cor.gc",
            "Difference profile\ndonor minus recipient, on aligned windows",
            "Circular binary segmentation\nchange points along the profile",
            "Putative introgressions\nsegments whose difference passes the threshold"
        ),
        kind = c("input", rep("step", 5), "result"),
        y = 7:1,
        stringsAsFactors = FALSE
    )
    nodes$x <- 1
    edges <- data.frame(from = 1:6, to = 2:7)
    g <- tidygraph::tbl_graph(nodes = nodes, edges = edges, node_key = "name")
    lay <- ggraph::create_layout(g, layout = "manual", x = nodes$x, y = nodes$y)

    ggraph::ggraph(lay) +
        ggraph::geom_edge_link(
            arrow = arrow(length = unit(2.6, "mm"), type = "closed"),
            end_cap = ggraph::circle(2.5, "mm"),
            colour = "grey55", edge_width = 0.5
        ) +
        ggraph::geom_node_label(
            aes(label = .data$label, fill = .data$kind),
            colour = "grey15", size = 2.9, lineheight = 1.05,
            label.r = unit(2, "mm"), label.padding = unit(2.4, "mm")
        ) +
        scale_fill_manual(
            name = NULL,
            values = c("input" = "#dbe9f6", "step" = "grey96", "result" = "#fdecea")
        ) +
        scale_x_continuous(limits = c(0.4, 1.6)) +
        scale_y_continuous(limits = c(0.3, 7.7)) +
        labs(
            title = "What MaMaMIA does to the sequencing data",
            caption = "Each step is a function in the package; the last one is the call to inspect."
        ) +
        theme_void(base_size = 11) +
        theme(
            plot.title = ggtext::element_markdown(face = "bold", size = 13, hjust = 0.5),
            plot.caption = ggtext::element_markdown(colour = "grey40", hjust = 0.5),
            plot.margin = margin(8, 8, 8, 8),
            legend.position = "none"
        )
}

# --------------------------------------------------------- 3. circos genome --

fig_circos <- function() {
    regions <- data.frame(
        chr = c(don_names, rec_names),
        start = 0,
        end = c(
            chrs$len[match(pairs$chr_id.don, chrs$chr_id)],
            chrs$len[match(pairs$chr_id.rec, chrs$chr_id)]
        )
    )
    diff_all <- do.call(rbind, lapply(seq_along(don_names), function(i) {
        out <- pair_out(i)
        rbind(
            data.frame(
                chr = don_names[i], start = out$start.don,
                end = out$stop.don, value = out$diff
            ),
            data.frame(
                chr = rec_names[i], start = out$start.rec,
                end = out$stop.rec, value = out$diff
            )
        )
    }))
    lim <- stats::quantile(abs(diff_all$value), 0.99, na.rm = TRUE)
    diff_all$value <- pmax(pmin(diff_all$value, lim), -lim)
    col_fun <- colorRamp2(c(-lim, 0, lim), c("#2166ac", "#f7f7f7", "#b2182b"))
    called <- isa$segments[isa$segments$putative_introgression, , drop = FALSE]

    f <- file.path(OUT, "03-circos.png")
    agg_png(f, width = 2000, height = 1800, res = 200, background = "white")
    on.exit(invisible(dev.off()), add = TRUE)
    par(mar = c(1, 1, 1, 1))
    circos.par(
        start.degree = 180, gap.degree = 2, points.overflow.warning = FALSE,
        track.margin = c(0.004, 0.004), cell.padding = c(0, 0, 0, 0)
    )
    circos.genomicInitialize(regions, plotType = NULL)
    circos.track(
        ylim = c(0, 1), bg.border = NA, track.height = mm_h(14),
        panel.fun = function(x, y) {
            s <- CELL_META$sector.index
            circos.rect(CELL_META$xlim[1], 0.15, CELL_META$xlim[2], 0.85,
                col = if (s %in% don_names) "#737373" else "#c9c9c9",
                border = "grey30", lwd = 0.4
            )
            circos.text(mean(CELL_META$xlim), 0.5, s,
                facing = "bending.inside", niceFacing = TRUE,
                cex = 0.7, col = "white"
            )
        }
    )
    circos.genomicTrack(
        diff_all, ylim = c(-lim, lim), track.height = mm_h(26),
        bg.border = "grey85",
        panel.fun = function(region, value, ...) {
            circos.lines(CELL_META$xlim, c(0, 0), col = "grey55", lwd = 0.5)
            circos.genomicRect(region, value,
                ytop.column = 1, ybottom = 0,
                col = col_fun(value[[1]]), border = NA
            )
        }
    )
    for (i in seq_len(nrow(called))) {
        s <- called[i, ]
        out <- isa$out[isa$out$chr_id.don == s$chr_id.don &
            isa$out$chr_id.rec == s$chr_id.rec, , drop = FALSE]
        circos.link(
            chrs$chr_name[match(s$chr_id.don, chrs$chr_id)],
            c(
                out$start.don[match(s$loc.start, out$iid)],
                out$stop.don[match(s$loc.end, out$iid)]
            ),
            chrs$chr_name[match(s$chr_id.rec, chrs$chr_id)],
            c(
                out$start.rec[match(s$loc.start, out$iid)],
                out$stop.rec[match(s$loc.end, out$iid)]
            ),
            col = "#cb181d33", border = "#cb181d"
        )
    }
    circos.clear()
    grid::grid.text(
        "At/A chromosomes: coverage difference and the segments called",
        y = 0.975, gp = grid::gpar(fontsize = 13, fontface = "bold")
    )
    grid::grid.text(
        "outer ring: donor (dark) and recipient (light) chromosomes    |    inner ring: donor - recipient coverage    |    ribbons: putative introgressions",
        y = 0.945, gp = grid::gpar(fontsize = 8, col = "grey35")
    )
    cat("wrote", f, "\n")
}

# ------------------------------------------------------------- 4. ideogram ---


# ---------------------------------------------------------------- 5. GC bias -

fig_gc_bias <- function() {
    subs <- unique(as.character(cd$data$subgenome))
    subs <- subs[order(subs)]
    d <- cd$data[cd$data$ideal, , drop = FALSE]
    d$subgenome <- as.character(d$subgenome)

    binned <- function(dd) {
        br <- seq(min(dd$gc), max(dd$gc), length.out = 13)
        b <- cut(dd$gc, breaks = br, include.lowest = TRUE)
        data.frame(
            gc = as.numeric(tapply(dd$gc, b, mean)),
            value = as.numeric(tapply(dd$cov, b, mean)),
            se = as.numeric(tapply(dd$cov, b, function(v) {
                stats::sd(v) / sqrt(length(v))
            })),
            n = as.integer(tapply(dd$cov, b, length)),
            subgenome = dd$subgenome[1],
            stringsAsFactors = FALSE
        )
    }
    pts <- do.call(rbind, lapply(subs, function(s) binned(d[d$subgenome == s, ])))
    pts$subgenome <- factor(pts$subgenome, levels = subs)

    grid <- expand.grid(
        gc = seq(min(d$gc), max(d$gc), length.out = 200L),
        subgenome = subs, stringsAsFactors = FALSE
    )
    grid$value <- MaMaMIA:::eta_zinb_response(cd$fit, grid)
    grid$subgenome <- factor(grid$subgenome, levels = subs)
    grid$donor <- grid$subgenome %in% c("At", "G")

    ggplot(pts, aes(x = .data$gc, y = .data$value)) +
        geom_errorbar(
            aes(ymin = .data$value - .data$se, ymax = .data$value + .data$se),
            width = 0.0012, colour = "grey45", linewidth = 0.35
        ) +
        geom_point(aes(size = .data$n), colour = "grey15", fill = "white", shape = 21) +
        geom_line(
            data = grid, aes(x = .data$gc, y = .data$value),
            colour = PAL$don, linewidth = 0.9
        ) +
        scale_size_continuous(
            name = "windows per bin", range = c(0.8, 3.4),
            breaks = c(50, 500, 1500)
        ) +
        facet_wrap(~subgenome, scales = "free_y", nrow = 1) +
        labs(
            title = "The GC response the model removes",
            subtitle = "Binned ideal-window coverage versus the fitted ZINB mean; all subgenomes peak near 45-46% GC.",
            caption = paste(
                "correctReadCounts() divides each window by this fitted response,",
                "putting equal true coverage on one scale."
            ),
            x = "GC content of the window", y = "Coverage"
        ) +
        theme_illus() +
        theme(
            legend.position = "right",
            strip.text = element_text(face = "bold"),
            plot.subtitle = ggtext::element_markdown(size = 9.2, lineheight = 1.05)
        )
}

# -------------------------------------------------------- 6. profile -> call -

fig_profile <- function(i = 1L) {
    out <- pair_out(i)
    seg <- pair_seg(i)
    seg$value <- seg[[isa$meta$metrics]]
    seg$x0 <- out$start.don[match(seg$loc.start, out$iid)] / 1e6
    seg$x1 <- out$stop.don[match(seg$loc.end, out$iid)] / 1e6
    seg$lwd <- ifelse(seg$putative_introgression, 2.4, 0.4)
    pts <- data.frame(x = out$start.don / 1e6, diff = out$diff)

    ggplot(pts, aes(x = .data$x, y = .data$diff)) +
        geom_hline(yintercept = 0, colour = "grey55", linetype = "dashed") +
        geom_point(size = 0.25, colour = "grey35", alpha = 0.7) +
        geom_segment(
            data = seg,
            aes(
                x = .data$x0, xend = .data$x1,
                y = .data$value, yend = .data$value,
                colour = .data$putative_introgression, linewidth = .data$lwd
            )
        ) +
        scale_linewidth_identity() +
        scale_colour_manual(
            name = NULL,
            values = c("FALSE" = "grey60", "TRUE" = PAL$don),
            labels = c("CBS segment", "putative introgression")
        ) +
        labs(
            title = paste0(
                "Segmentation of one pair: *", don_names[i], "* / ", rec_names[i]
            ),
            subtitle = paste0(
                "Bars are the ", isa$meta$metrics,
                " of each CBS segment; the red one passes the threshold"
            ),
            x = "Donor position (Mb)", y = "Donor - recipient coverage"
        ) +
        theme_illus()
}

# ------------------------------------------------------------ 7. orientation -

fig_orientation <- function(i = 7L) {
    chr_id <- pairs$chr_id.don[i]
    rev_isa <- MaMaMIA::reverseWindows(isa, chr_ids = chr_id)
    collect <- function(obj, label) {
        out <- obj$out[obj$out$chr_id.don == pairs$chr_id.don[i] &
            obj$out$chr_id.rec == pairs$chr_id.rec[i], , drop = FALSE]
        seg <- obj$segments[obj$segments$chr_id.don == pairs$chr_id.don[i] &
            obj$segments$chr_id.rec == pairs$chr_id.rec[i], , drop = FALSE]
        seg$value <- seg[[obj$meta$metrics]]
        list(
            pts = data.frame(iid = out$iid, diff = out$diff, orientation = label),
            seg = data.frame(
                loc.start = seg$loc.start, loc.end = seg$loc.end,
                value = seg$value, called = seg$putative_introgression,
                orientation = label
            )
        )
    }
    a <- collect(isa, "as assembled")
    b <- collect(rev_isa, "donor chromosome reversed")
    a$seg$lwd <- ifelse(a$seg$called, 2.4, 0.4)
    b$seg$lwd <- ifelse(b$seg$called, 2.4, 0.4)
    pts <- rbind(a$pts, b$pts)
    seg <- rbind(a$seg, b$seg)
    lv <- c("as assembled", "donor chromosome reversed")
    pts$orientation <- factor(pts$orientation, levels = lv)
    seg$orientation <- factor(seg$orientation, levels = lv)

    ggplot(pts, aes(x = .data$iid, y = .data$diff)) +
        geom_hline(yintercept = 0, colour = "grey55", linetype = "dashed") +
        geom_point(size = 0.25, colour = "grey35", alpha = 0.7) +
        geom_segment(
            data = seg,
            aes(
                x = .data$loc.start, xend = .data$loc.end,
                y = .data$value, yend = .data$value,
                colour = .data$called, linewidth = .data$lwd
            )
        ) +
        scale_linewidth_identity() +
        scale_colour_manual(
            name = NULL,
            values = c("FALSE" = "grey60", "TRUE" = PAL$don),
            labels = c("CBS segment", "putative introgression")
        ) +
        facet_wrap(~orientation, ncol = 1) +
        labs(
            title = paste0(
                "Assembly orientation moves the signal: *",
                don_names[i], "* / ", rec_names[i]
            ),
            subtitle = paste0(
                "The donor coverage vector is reversed before segmentation; calls change from ",
                sum(a$seg$called), "/", nrow(a$seg), " to ",
                sum(b$seg$called), "/", nrow(b$seg), " segments."
            ),
            x = "Window along the chromosome", y = "Donor - recipient coverage"
        ) +
        theme_illus() +
        theme(plot.subtitle = ggtext::element_markdown(size = 9.2, lineheight = 1.05))
}


# ------------------------------------------------- 9. arrow map (gggenes) ----

fig_arrowmap <- function() {
    rows <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
        seg <- pair_seg(i)
        called <- seg[seg$putative_introgression, , drop = FALSE]
        out <- pair_out(i)
        if (nrow(called) > 0L) {
            called$s <- out$start.don[match(called$loc.start, out$iid)]
            called$e <- out$stop.don[match(called$loc.end, out$iid)]
        }
        side <- function(nm, len) {
            body <- data.frame(
                chr = nm, start = 0, end = len, forward = TRUE,
                type = "chromosome", stringsAsFactors = FALSE
            )
            if (nrow(called) == 0L) {
                return(body)
            }
            rbind(body, data.frame(
                chr = nm, start = called$s, end = called$e, forward = TRUE,
                type = "putative introgression", stringsAsFactors = FALSE
            ))
        }
        rbind(
            side(don_names[i], chrs$len[match(pairs$chr_id.don[i], chrs$chr_id)]),
            side(rec_names[i], chrs$len[match(pairs$chr_id.rec[i], chrs$chr_id)])
        )
    }))
    rows$chr <- factor(rows$chr, levels = unique(rows$chr))
    rows$type <- factor(rows$type,
        levels = c("chromosome", "putative introgression")
    )

    ggplot(rows, aes(
        xmin = .data$start / 1e6, xmax = .data$end / 1e6,
        y = .data$chr, fill = .data$type, forward = .data$forward
    )) +
        geom_gene_arrow(
            arrowhead_height = unit(2.2, "mm"), arrowhead_width = unit(2.2, "mm")
        ) +
        scale_fill_manual(
            name = NULL,
            values = c(
                "chromosome" = "grey80",
                "putative introgression" = PAL$don
            )
        ) +
        theme_genes() +
        labs(
            title = "Introgression map: each chromosome as a gene-arrow track",
            subtitle = paste0(
                "The same 14 chromosomes as in the other views; red arrows are the ",
                "segments passing the threshold"
            ),
            x = "Position (Mb)", y = NULL
        ) +
        theme(
            plot.title = ggtext::element_markdown(face = "bold", size = 12),
            plot.subtitle = ggtext::element_markdown(colour = "grey30", size = 9),
            legend.position = "bottom"
        )
}

# ------------------------------------------- 10. distributions (ggdist) -----

fig_distributions <- function() {
    diffs <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
        out <- pair_out(i)
        data.frame(
            pair = paste0(don_names[i], " / ", rec_names[i]),
            diff = out$diff
        )
    }))
    ord <- names(sort(tapply(diffs$diff, diffs$pair, stats::median)))
    diffs$pair <- factor(diffs$pair, levels = ord)

    ggplot(diffs, aes(x = .data$diff, y = .data$pair)) +
        ggdist::stat_halfeye(
            .width = c(0.5, 0.95), fill = PAL$rec, colour = "grey30",
            slab_alpha = 0.75, point_colour = "grey15", point_size = 1.1,
            normalize = "groups"
        ) +
        geom_vline(xintercept = 0, colour = "grey40", linetype = "dashed") +
        labs(
            title = "Where the donor-recipient difference actually sits",
            subtitle = paste0(
                "Distribution of the per-window difference for each pair, ordered by ",
                "median; the point is the median with 50% and 95% intervals"
            ),
            x = "Donor - recipient coverage", y = NULL
        ) +
        theme_illus()
}

# -------------------------------------------------------------------- main ---

save_fig(fig_mechanism(), "01-mechanism.png", 9, 6)
save_fig(fig_workflow(), "02-workflow.png", 8.5, 7)
fig_circos()
save_fig(fig_gc_bias(), "04-gc-bias.png", 9, 6.5)
save_fig(fig_profile(), "05-pair-profile.png", 9, 4)
save_fig(fig_orientation(), "06-orientation.png", 9, 6)
save_fig(fig_arrowmap(), "07-arrowmap.png", 9, 7)
save_fig(fig_distributions(), "08-distributions.png", 9, 5)

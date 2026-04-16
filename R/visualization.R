#' Overview of cellular networks with a bubble plot
#'
#' @param obj  A SCSRNet or SCSRNoNet object.
#' @param selected.populations A vector of cell population names
#' to consider in the plot. By default, all the populations are
#' considered.
#' @param genes.to.count A vector of gene names for counting or enrichment
#' analysis.
#' @param only.R.in.genes  A logical indicating whether \code{genes.to.count}
#' should be regarded as containing receptor gene names only.
#' @param use.proportions  A logical to choose between representing the
#' proportion of genes in \code{genes.to.count} or its enrichment.
#' @param low.color  The color to be used when no gene list is provided
#' or when the proportion of genes in \code{genes.to.count} is 0.
#' @param high.color  The color for maximum proportion or best enrichment
#' P-value.
#'
#' @details A matrix dot plot is generated to represent how much each pair of
#' cell populations interact including autocrine interactions when
#' available.
#'
#' By default, the plot only reports the number of interactions, but it is
#' possible to provide a list of genes of interest (\code{genes.to.count}
#' parameter). One option in this
#' case is to provide receptors involved in specific signaling whose
#' abundance is to be illustrated on top of the number of interactions. For
#' this, \code{only.R.in.genes} must be set to TRUE, the default. If not,
#' then only interactions with both the ligand and the receptor in
#' \code{genes.to.count} are considered. This may enable more specific
#' counting. Lastly, it is possible to choose between color-coding the
#' proportion of interactions with receptor or both receptor and ligand in
#' \code{genes.to.count}, or to perform an enrichment analysis. In the
#' last case, the color-code is based on -log10(P-values).
#'
#' @return A bubble plot is displayed reflecting population pairwise
#' interactions.
#'
#' @export
#'
#' @import ggplot2
#' @importFrom BulkSignalR LRinter
#' @examples
#' print("cellNetBubblePlot")
#' if (FALSE) {
#'     cellNetBubblePlot(scsrcn)
#' }
cellNetBubblePlot <- function(obj, selected.populations = NULL,
    genes.to.count = NULL, only.R.in.genes = TRUE,
    use.proportions = FALSE,
    low.color = "gray25", high.color = "firebrick1") {
    
    from <- to <- col.fact <- NULL

    if (!is(obj, "SCSRNoNet") && !is(obj, "SCSRNet")) {
        stop("obj must be either of SCSRNoNet or SCSRNet class")
    }
    if (!is.null(selected.populations) & 
        !all(selected.populations %in% populations(obj))) {
        stop("selected.populations must all be in obj populations")
    }

    cell.net <- is(obj, "SCSRNet")
    if (is.null(selected.populations)) {
        selected.populations <- unique(populations(obj))
    }
    auto <- autocrines(obj)
    para <- paracrines(obj)

    # background gene list for hypergeometric test
    bg.genes <- NULL
    if (!is.null(genes.to.count)) {
        if (length(auto) > 0) {
            for (pop in selected.populations) {
                if (pop %in% names(auto)) {
                    if (cell.net) {
                        tab <- BulkSignalR::LRinter(auto[[pop]])[,
                        c("L", "R")]
                    } else {
                        tab <- auto[[pop]][, c("L", "R")]
                    }
                    if (only.R.in.genes) {
                        bg.genes <- c(bg.genes, tab$R)
                    } else {
                        bg.genes <- c(bg.genes, tab$L, tab$R)
                    }
                }
            }
        }
        if (length(para) > 0) {
            for (source.pop in selected.populations) {
                for (target.pop in selected.populations) {
                    if (source.pop != target.pop) {
                        p <- paste0(source.pop, "_vs_", target.pop)
                        if (p %in% names(para)) {
                            if (cell.net) {
                                tab <- BulkSignalR::LRinter(para[[p]])[,
                                c("L", "R")]
                            } else {
                                tab <- para[[p]][, c("L", "R")]
                            }
                            if (only.R.in.genes) {
                                bg.genes <- c(bg.genes, tab$R)
                            } else {
                                bg.genes <- c(bg.genes, tab$L, tab$R)
                            }
                        }
                    }
                }
            }
        }
        bg.genes <- unique(bg.genes)
    }

    # count interactions and prepare data.frame for ggplot
    comm <- NULL
    if (length(auto) > 0) {
        for (pop in selected.populations) {
            if (pop %in% names(auto)) {
                if (cell.net) {
                    tab <- unique(BulkSignalR::LRinter(auto[[pop]])[,
                    c("L", "R")])
                } else {
                    tab <- unique(auto[[pop]][, c("L", "R")])
                }
                n <- nrow(tab)
                if (is.null(genes.to.count)) {
                    n.in.genes <- 0
                    prop <- 0
                    pval <- 1
                } else {
                    n.in.genes <- ifelse(only.R.in.genes,
                        sum(tab$R %in% genes.to.count),
                        sum(tab$L %in% genes.to.count & 
                            tab$R %in% genes.to.count)
                    )
                    prop <- n.in.genes / n
                    k <- ifelse(only.R.in.genes,
                        length(unique(tab$R)),
                        nrow(tab)
                    )
                    pval <- stats::phyper(
                        q = n.in.genes, 
                        m = length(intersect(genes.to.count, bg.genes)),
                        n = length(setdiff(bg.genes, genes.to.count)),
                        k = k, lower.tail = FALSE
                    )
                }
                comm <- rbind(comm,
                    data.frame(from = pop, 
                        to = pop, 
                        n = n, 
                        prop = prop, 
                        pval = pval))
            }
        }
    }
    if (length(para) > 0) {
        for (source.pop in selected.populations) {
            for (target.pop in selected.populations) {
                if (source.pop != target.pop) {
                    p <- paste0(source.pop, "_vs_", target.pop)
                    if (p %in% names(para)) {
                        if (cell.net) {
                            tab <- unique(BulkSignalR::LRinter(para[[p]])[,
                                c("L", "R")])
                        } else {
                            tab <- unique(para[[p]][, c("L", "R")])
                        }
                        n <- nrow(tab)
                        if (is.null(genes.to.count)) {
                            n.in.genes <- 0
                            prop <- 0
                            pval <- 1
                        } else {
                            n.in.genes <- ifelse(only.R.in.genes,
                                sum(tab$R %in% genes.to.count),
                                sum(tab$L %in% genes.to.count & 
                                    tab$R %in% genes.to.count)
                            )
                            prop <- n.in.genes / n
                            k <- ifelse(only.R.in.genes,
                                length(unique(tab$R)),
                                nrow(tab)
                            )
                            pval <- stats::phyper(
                                q = n.in.genes, 
                                m = length(intersect(genes.to.count, bg.genes)),
                                n = length(setdiff(bg.genes, genes.to.count)), 
                                k = k,
                                lower.tail = FALSE
                            )
                        }
                        comm <- rbind(comm, data.frame(from = source.pop,
                            to = target.pop,
                            n = n,
                            prop = prop,
                            pval = pval))
                    }
                }
            }
        }
    }
    comm$from <- as.factor(comm$from)
    comm$to <- as.factor(comm$to)

    # produce the plot
    if (is.null(genes.to.count)) {
        g <- ggplot2::ggplot(comm, ggplot2::aes(x = from, y = to)) +
            ggplot2::geom_point(ggplot2::aes(size = n)) +
            ggplot2::scale_size(name = "# inter") +
            ggplot2::theme_bw() +
            ggplot2::theme(axis.text.x = 
                ggplot2::element_text(angle = 45, hjust = 1))
    } else {
        if (use.proportions) {
            comm$col.fact <- comm$prop
        } else {
            comm$col.fact <- -log10(comm$pval)
        }
        g <- ggplot2::ggplot(comm, ggplot2::aes(x = from, y = to)) +
            ggplot2::geom_point(ggplot2::aes(color = col.fact, size = n)) +
            ggplot2::scale_size(name = "# inter") +
            ggplot2::scale_color_gradient(low = low.color, high = high.color) +
            ggplot2::theme_bw() +
            ggplot2::theme(axis.text.x = 
                ggplot2::element_text(angle = 45, hjust = 1)) +
            ggplot2::labs(color = 
                ifelse(use.proportions, "Proportion", "-log10(P)"))
    }

    plot(g)

} # cellNetBubblePlot


#' Heatmap overview of cellular networks
#'
#' @param obj  A SCSRNet or SCSRNoNet object.
#' @param selected.populations A vector of cell population names
#' to consider in the plot. By default, all the populations are
#' considered.
#' @param genes.to.count A vector of gene names for counting or enrichment
#' analysis.
#' @param only.R.in.genes  A logical indicating whether \code{genes.to.count}
#' should be regarded as containing receptor gene names only.
#' @param use.proportions  A logical to choose between representing the
#' proportion of genes in \code{genes.to.count} or its enrichment.
#' @param low.color  The color to be used for the lowest value.
#' @param high.color  The color to be used for the highest value.
#' @param thres A higher threshold imposed to the values.
#' @param col.fun A function that returns a color based on the values in the
#' matrix. When no such function is provided, a function is generated from
#' \code{low.color} and \code{high.color} with a linear gradient.
#' @param row.bar.color The color of the barplot counting the total number of
#' received LR interactions.
#' @param col.bar.color The color of the barplot counting the total number of
#' emitted LR interactions.
#' @param pop.font.size Font size for cell population names.
#' @param title.font.size Font size for row and column titles.
#' @param legend.font.size Font size for the legends.
#' @param bar.plot.height Heigh (or width) or the two barplots.
#' @param bar.num A logical indicating whether numbers should be written on top
#' of the barplots.
#' @param bar.num.font.size Font size of the totals on the barplots.
#' @param original.order A logical indicating whether matrix rows and columns
#' should be left in their original order. If not (the default), they are
#' reordered by computing two dendrograms that are not displayed.
#' @param rect.border.width White border width around the colored
#' rectangles of the heatmap.
#'
#' @details A heatmap is generated to represent how much each pair of
#' cell populations interact including autocrine interactions when
#' available.
#'
#' By default, the plot only reports the number of interactions, but it is
#' possible to provide a list of genes of interest (\code{genes.to.count}
#' parameter). One option in this
#' case is to provide receptors involved in specific signaling whose
#' abundance is to be illustrated on top of the number of interactions. For
#' this, \code{only.R.in.genes} must be set to TRUE, the default. If not,
#' then only interactions with both the ligand and the receptor in
#' \code{genes.to.count} are considered. This may enable more specific
#' counting. Lastly, it is possible to choose between color-coding the
#' proportion of interactions with receptor or both receptor and ligand in
#' \code{genes.to.count}, or to perform an enrichment analysis. In the
#' last case, the color-code is based on -log10(P-values).
#'
#' @return A heatmap is displayed reflecting population pairwise interactions.
#'
#' @export
#'
#' @import ComplexHeatmap
#' @import circlize
#' @importFrom grid gpar
#' @importFrom BulkSignalR LRinter
#' @examples
#' print("cellNetHeatmap")
#' if (FALSE) {
#'     cellNetHeatmap(scsrcn)
#' }
cellNetHeatmap <- function(obj, selected.populations = NULL,
                           genes.to.count = NULL, only.R.in.genes = TRUE,
                           use.proportions = FALSE,
                           low.color = "white", high.color = "royalblue3",
                           thres = NULL, col.fun = NULL,
                           row.bar.color = "darkcyan",
                           col.bar.color = "seagreen",
                           pop.font.size = 10,
                           title.font.size = 12, legend.font.size = 10,
                           bar.plot.height = NULL, rect.border.width = 4,
                           bar.num = FALSE, bar.num.font.size = 8,
                           original.order = FALSE){
    
    if (!is(obj, "SCSRNoNet") && !is(obj, "SCSRNet")) {
        stop("obj must be either of SCSRNoNet or SCSRNet class")
    }
    if (!is.null(selected.populations) & 
        !all(selected.populations %in% populations(obj))) {
        stop("selected.populations must all be in obj populations")
    }
    
    cell.net <- is(obj, "SCSRNet")
    if (is.null(selected.populations)) {
        selected.populations <- unique(populations(obj))
    }
    auto <- autocrines(obj)
    para <- paracrines(obj)
    
    # background gene list for hypergeometric test
    bg.genes <- NULL
    if (!is.null(genes.to.count)) {
        if (length(auto) > 0) {
            for (pop in selected.populations) {
                if (pop %in% names(auto)) {
                    if (cell.net) {
                        tab <- BulkSignalR::LRinter(auto[[pop]])[,
                                                                 c("L", "R")]
                    } else {
                        tab <- auto[[pop]][, c("L", "R")]
                    }
                    if (only.R.in.genes) {
                        bg.genes <- c(bg.genes, tab$R)
                    } else {
                        bg.genes <- c(bg.genes, tab$L, tab$R)
                    }
                }
            }
        }
        if (length(para) > 0) {
            for (source.pop in selected.populations) {
                for (target.pop in selected.populations) {
                    if (source.pop != target.pop) {
                        p <- paste0(source.pop, "_vs_", target.pop)
                        if (p %in% names(para)) {
                            if (cell.net) {
                                tab <- BulkSignalR::LRinter(para[[p]])[,
                                                                       c("L", "R")]
                            } else {
                                tab <- para[[p]][, c("L", "R")]
                            }
                            if (only.R.in.genes) {
                                bg.genes <- c(bg.genes, tab$R)
                            } else {
                                bg.genes <- c(bg.genes, tab$L, tab$R)
                            }
                        }
                    }
                }
            }
        }
        bg.genes <- unique(bg.genes)
    }
    
    # count interactions and prepare a matrix for the heatmap
    n.pop <- length(selected.populations)
    comm.n <- matrix(0, nrow=n.pop, ncol=n.pop,
                     dimnames=list(selected.populations, selected.populations)
    )
    comm.prop <- comm.n
    comm.pval <- matrix(1, nrow=n.pop, ncol=n.pop,
                    dimnames=list(selected.populations, selected.populations)
    )
    if (length(auto) > 0) {
        for (pop in selected.populations) {
            if (pop %in% names(auto)) {
                if (cell.net) {
                    tab <- unique(BulkSignalR::LRinter(auto[[pop]])[,
                                                                    c("L", "R")])
                } else {
                    tab <- unique(auto[[pop]][, c("L", "R")])
                }
                n <- nrow(tab)
                if (is.null(genes.to.count)) {
                    n.in.genes <- 0
                    prop <- 0
                    pval <- 1
                } else {
                    n.in.genes <- ifelse(only.R.in.genes,
                                         sum(tab$R %in% genes.to.count),
                                         sum(tab$L %in% genes.to.count & 
                                                 tab$R %in% genes.to.count)
                    )
                    prop <- n.in.genes / n
                    k <- ifelse(only.R.in.genes,
                                length(unique(tab$R)),
                                nrow(tab)
                    )
                    pval <- stats::phyper(
                        q = n.in.genes, 
                        m = length(intersect(genes.to.count, bg.genes)),
                        n = length(setdiff(bg.genes, genes.to.count)),
                        k = k, lower.tail = FALSE
                    )
                }
                comm.n[pop, pop] <- n
                comm.prop[pop, pop] <- prop
                comm.pval[pop, pop] <- pval
            }
        }
    }
    if (length(para) > 0) {
        for (source.pop in selected.populations) {
            for (target.pop in selected.populations) {
                if (source.pop != target.pop) {
                    p <- paste0(source.pop, "_vs_", target.pop)
                    if (p %in% names(para)) {
                        if (cell.net) {
                            tab <- unique(BulkSignalR::LRinter(para[[p]])[,
                                                                          c("L", "R")])
                        } else {
                            tab <- unique(para[[p]][, c("L", "R")])
                        }
                        n <- nrow(tab)
                        if (is.null(genes.to.count)) {
                            n.in.genes <- 0
                            prop <- 0
                            pval <- 1
                        } else {
                            n.in.genes <- ifelse(only.R.in.genes,
                                                 sum(tab$R %in% genes.to.count),
                                                 sum(tab$L %in% genes.to.count & 
                                                         tab$R %in% genes.to.count)
                            )
                            prop <- n.in.genes / n
                            k <- ifelse(only.R.in.genes,
                                        length(unique(tab$R)),
                                        nrow(tab)
                            )
                            pval <- stats::phyper(
                                q = n.in.genes, 
                                m = length(intersect(genes.to.count, bg.genes)),
                                n = length(setdiff(bg.genes, genes.to.count)), 
                                k = k,
                                lower.tail = FALSE
                            )
                        }
                        comm.n[target.pop, source.pop] <- n
                        comm.prop[target.pop, source.pop] <- prop
                        comm.pval[target.pop, source.pop] <- pval
                    }
                }
            }
        }
    }

    # produce the plot --------------------
    
    # type of plot and right data into m
    if (is.null(genes.to.count)) {
        m <- comm.n
        plot.type <- "num"
    } else {
        if (use.proportions) {
            m <- comm.prop
            plot.type <- "prop"
        } else {
            m <- -log10(comm.pval)
            m[!is.finite(m)] <- 0
            plot.type <- "pval"
        }
    }
        
    # data preparation
    if (original.order){
        dend.row = FALSE
        dend.col = FALSE
    }
    else{
        dist.row <- stats::dist(m)
        dend.row <- stats::as.dendrogram(stats::hclust(dist.row,
                                                       method="ward.D"))
        dist.col <- stats::dist(t(m))
        dend.col <- stats::as.dendrogram(stats::hclust(dist.col,
                                                       method="ward.D"))
    }
    max0 <- max(m)
    if (!is.null(thres)){
        m[m > thres] <- thres
    }
    if (is.null(col.fun)){
        col.fun <- circlize::colorRamp2(breaks=c(0, max(m), max0),
                              c(low.color, high.color, high.color))
    }
    
    if (plot.type == "num"){
        # side annotations
        top.annot <- ComplexHeatmap::HeatmapAnnotation(
            emissions = ComplexHeatmap::anno_barplot(colSums(m), border=FALSE,
                                 gp=gpar(fill=col.bar.color, col=col.bar.color),
                                 add_numbers=bar.num,
                                 numbers_gp=gpar(fontsize=bar.num.font.size)),
            height = bar.plot.height,
            show_annotation_name = FALSE
        )
        right.annot <- ComplexHeatmap::rowAnnotation(
            receptions = ComplexHeatmap::anno_barplot(rowSums(m), border=FALSE,
                                 gp=gpar(fill=row.bar.color, col=row.bar.color),
                                 which="row", add_numbers=bar.num,
                                 numbers_gp=gpar(fontsize=bar.num.font.size)),
            width = bar.plot.height,
            show_annotation_name = FALSE
        )
    
        # plot
        ComplexHeatmap::Heatmap(m, col=col.fun, cluster_rows=dend.row,
            cluster_columns=dend.col, show_row_dend=FALSE,
            show_column_dend=FALSE, row_names_side="left",
            rect_gp=gpar(col="white", lwd=rect.border.width),
            column_title="From", column_title_gp=gpar(fontsize=title.font.size,
                                                  fontface="bold"),
            row_title="To", row_title_gp=gpar(fontsize=title.font.size,
                                              fontface="bold"),
            row_names_gp=gpar(fontsize=pop.font.size),
            column_names_gp=gpar(fontsize=pop.font.size),
            top_annotation=top.annot, right_annotation=right.annot,
            heatmap_legend_param=list(title="# LRIs", fontsize=legend.font.size)
        )
    }
    else{
        ComplexHeatmap::Heatmap(m, col=col.fun, cluster_rows=dend.row,
            cluster_columns=dend.col, show_row_dend=FALSE,
            show_column_dend=FALSE, row_names_side="left",
            rect_gp=gpar(col="white", lwd=rect.border.width),
            column_title="From", column_title_gp=gpar(fontsize=title.font.size,
                                                      fontface="bold"),
            row_title="To", row_title_gp=gpar(fontsize=title.font.size,
                                              fontface="bold"),
            row_names_gp=gpar(fontsize=pop.font.size),
            column_names_gp=gpar(fontsize=pop.font.size),
            heatmap_legend_param=list(title=ifelse(use.proportions,
                                                   "Proportion",
                                                   "-log10(P-value)"),
                                      fontsize=legend.font.size)
        )
    }
    
} # cellNetHeatmap

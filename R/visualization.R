#' Overview of cellular networks
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
#' @export
#'
#' @return A BSRDataModel object containing the smoothed ncounts.
#'
#' @import ggplot2
#' @importFrom BulkSignalR LRinter
#' @examples
#' print("cellNetBubblePlot")
#' if (FALSE) {
#'     cellBubblePlot(scsrcn)
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

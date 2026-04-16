#' SingleCellSignalR no network object
#'
#' An S4 class to represent data and naive, network-free inferences
#' of ligand-receptor interactions. Both
#' autocrine and paracrine interactions can be inferred and explored.
#' The absence of receptor network downstream exploration makes computations
#' fast, but the quality of the inferences is worse than those
#' obtained using the SCSRNet class.
#'
#' @slot bsrdm.comp   A BulkSignalR BSRDataModelComp object containing
#' the expression matrix data as well as the comparisons between cell
#' populations.
#' @slot populations  A vector defining the cell population to which
#' each individual cell belongs to. The length of this vector must
#' be equal to the number of columns of the expression matrix.
#' @slot paracrines   A list of data.frames
#' containing all the paracrine ligand-receptor interactions between
#' all the possible pairs of cell populations.
#' @slot autocrines   A list of data.frames
#' containing all the autocrine ligand-receptor interactions for
#' each cell population.
#' @slot inf.param    Inference parameters.
#'
#' @details  This class is a container for all the data and inferences
#' related to a given single-cell project.
#' @export
#' @examples
#' print("Create SCSRNet object:")
#' data(example_dataset, package = "SingleCellSignalR")
#' mat <- log1p(data.matrix(example_dataset[, -1])) / log(2)
#' rownames(mat) <- example_dataset[[1]]
#' rme <- rowMeans(mat)
#' mmat <- mat[rme > 0.05, ]
#' d <- dist(t(mmat))
#' h <- hclust(d, method = "ward.D")
#' pop <- paste0("pop_", cutree(h, 5))
#' scsrnn <- SCSRNoNet(mat,
#'     populations = pop, normalize = FALSE,
#'     method = "log-only", log.transformed = TRUE
#' )
#' @importFrom BulkSignalR coerce
#' @importClassesFrom BulkSignalR BSRDataModel BSRDataModelComp
#' @importFrom methods is as
setClass("SCSRNoNet",
    slots = c(
        bsrdm.comp = "BSRDataModelComp",
        populations = "character",
        paracrines = "list",
        autocrines = "list",
        inf.param = "list"
    ),
    prototype = list(
        bsrdm.comp = as(
            new("BSRDataModel"),
            "BSRDataModelComp"
        ),
        populations = "B-cells",
        paracrines = list(),
        autocrines = list(),
        inf.param = list()
    )
)

setValidity(
    "SCSRNoNet",
    function(object) {
        if (!is.character(object@populations)) {
            return("populations is not of character type")
        }
        if (!is(object@bsrdm.comp, "BSRDataModelComp")) {
            return("bsrdm.comp is not of class BSRDataModelComp")
        }
        if (!is.list(object@paracrines)) {
            return("paracrines is not a list")
        }
        if (!is.list(object@autocrines)) {
            return("autocrines is not a list")
        }
        if (!is.list(object@inf.param)) {
            return("inf.param is not a list")
        }

        TRUE
    }
)

setMethod(
    "show", "SCSRNoNet",
    function(object) {
        print(object@bsrdm.comp)
        utils::str(object@inf.param)
        cat("Cell populations:\n")
        print(table(object@populations))
        cat(
            "Paracrine interactions:",
            paste(
                names(object@paracrines)[seq_len(min(
                    5,
                    length(names(object@paracrines))
                ))],
                collapse = ", "
            ),
            "\n"
        )
        cat(
            "Autocrine interactions:",
            paste(
                names(object@autocrines)[seq_len(min(
                    5,
                    length(names(object@autocrines))
                ))],
                collapse = ", "
            ),
            "\n"
        )
    }
)


# Constructor ========================================================

#' Prepare a SCSRNoNet object from expression data
#'
#' Take a matrix or data frame containing single-cell RNA sequencing
#' or proteomics data and return a SCSRNoNet
#' object ready for subsequent cellular network inference.
#' Normally, SCSRNoNet objects
#' are not instantiated directly, but through this function.
#'
#' @param counts     A table or matrix of read counts.
#' @param species    Data were obtained for this organism.
#' @param populations  A vector indicating to which cell population
#'   each individual cell belongs to.
#' @param normalize  A logical indicating whether \code{counts} should be
#'   normalized according to \code{method} or if it was normalized beforehand.
#' @param symbol.col The index of the column containing the gene symbols in case
#'   those are not the row names of \code{counts} already.
#' @param min.count  The minimum read count of a gene to be considered expressed
#'   in a sample.
#' @param prop       The minimum proportion of samples where a gene must be
#'   expressed higher than \code{min.count} to keep that gene.
#' @param method     The normalization method ('UQ' for upper quartile or 'TC'
#'   for total count). If \code{normalize==FALSE}, then method must be
#'   used to document the name of the normalization method applied by the user.
#' @param UQ.pc      Percentile for upper-quartile normalization, number
#' between 0 and 1 (in case the default 0.95 is not
#' appropriate).
#' @param log.transformed  A logical indicating whether expression data were
#'   already log2-transformed.
#' @param min.LR.found  The minimum number of ligands or receptors found in
#'   \code{count} row names after eliminating the rows containing too many
#'   zeros according to \code{min.count} and \code{prop}.
#' @param conversion.dict  Correspondence table of HUGO gene symbols
#' human/nonhuman. Not used unless the organism is different from human.
#'
#' @return A SCSRNoNet object with empty interactions.
#'
#' @details The \code{counts} matrix or table should be provided with expression
#'   levels of protein coding genes in each samples (column) and
#'   \code{rownames(counts)} set to HUGO official gene symbols. 
#'   For commodity, it is also possible to provide \code{counts} with the
#'   gene symbols stored in one of its columns. This column must be specified
#'   with \code{symbol.col}. In such a case, \code{prepareDataset} will extract
#'   this column and use it to set the row names. Because row names must be
#'   unique, \code{prepareDataset} will eliminate rows with duplicated gene
#'   symbols by keeping the rows with maximum average expression. Gene symbol
#'   duplication may occur in protein coding genes after genome alignment
#'   due to errors in genome feature annotation files (GTF/GFF), where a handful
#'   of deprecated gene annotations might remain, or
#'   some genes are not given their fully specific symbols. If your read count
#'   extraction pipeline does not take care of this phenomenon, the maximum mean
#'   expression selection strategy implemented here should solve this difficulty
#'   for the sake of inferring ligand-receptor interactions.
#'
#'   If \code{normalize} is \code{TRUE} then normalization is performed
#'   according to \code{method}. If those two simple methods are not satisfying,
#'   then it is possible to provide a pre-normalized matrix setting
#'   \code{normalize} to \code{FALSE}. In such a case, the parameter
#'   \code{method} must be used to document the name of the normalization
#'   algorithm used.
#'
#'   In case proteomic data are provided, \code{min.count} must be
#'   understood as its equivalent with respect to those data.
#'
#'
#' @export
#' @importFrom BulkSignalR coerce BSRDataModel
#' @importClassesFrom BulkSignalR BSRDataModel BSRDataModelComp
#' @importFrom methods as
#' @examples
#' print("SCSRNoNet")
#' data(example_dataset, package = "SingleCellSignalR")
#' mat <- log1p(data.matrix(example_dataset[, -1])) / log(2)
#' rownames(mat) <- example_dataset[[1]]
#' rme <- rowMeans(mat)
#' mmat <- mat[rme > 0.05, ]
#' d <- dist(t(mmat))
#' h <- hclust(d, method = "ward.D")
#' pop <- paste0("pop_", cutree(h, 5))
#' scsrnn <- SCSRNoNet(mat,
#'     normalize = FALSE, method = "log-only",
#'     log.transformed = TRUE, populations = pop
#' )
#'
SCSRNoNet <- function(counts, populations, normalize = TRUE,
    symbol.col = NULL, min.count = 1,
    prop = 0.1, method = c("UQ", "TC"), 
    log.transformed = FALSE, min.LR.found = 80,
    species = "hsapiens", conversion.dict = NULL,
    UQ.pc = 0.95) {
    if ((species != "hsapiens") && is.null(conversion.dict)) {
        stop("Non-human species but no conversion.dict provided")
    }
    if (length(populations) != ncol(counts)) {
        stop("The populations vector has a length != number of cells in counts")
    }

    if (normalize) {
        if (prop < 0 || prop > 1) {
            stop("prop must lie in [0;1]")
        }
        if (UQ.pc <= 0 || UQ.pc > 1) {
            stop("UQ.pc must lie in ]0;1]")
        }
        if (min.count < 0) {
            stop("min.count must be positive")
        }
        method <- match.arg(method)
    } else if (nchar(method) == 0) {
        stop(
            "In case of user-normalized counts, the name of the ",
            "normalization must be documented through the parameter 'method'"
        )
    }

    if (!is.null(symbol.col)) {
        if (!is.numeric(symbol.col)) {
            stop("symbol.col must be the index",
                " of the column containing the gene symbols")
        }

        # simple but desperately slow counts <-
        # aggregate(.~symbol,data=counts,FUN=max)

        # home-made but fast
        symbols <- as.character(counts[, symbol.col])
        d <- symbols[duplicated(symbols)]
        bad <- NULL
        for (s in d) {
            i <- which(symbols == s)
            t <- rowSums(counts[i, -symbol.col])
            bad <- c(bad, i[-which.max(t)])
        }

        # remove duplicates and the gene symbol column
        if (!is.null(bad)) {
            counts <- counts[-bad, -symbol.col]
            rownames(counts) <- symbols[-bad]
        } else {
            counts <- counts[, -symbol.col]
            rownames(counts) <- symbols
        }
    }

    if (is.null(rownames(counts)) || typeof(rownames(counts)) != "character") {
        stop("The read count matrix must ",
            " be provided with gene symbols as row names")
    }

    # as of now we ensure that counts is a matrix
    if (!is.matrix(counts)) {
        counts <- data.matrix(counts)
    }

    # avoid empty rows even if no normalization is performed here
    counts <- counts[rowSums(abs(counts)) > 0, ]

    if (normalize) {
        good.c <- rowSums(counts >= min.count) >= prop * ncol(counts)
        counts <- counts[good.c, ]
        if (method == "UQ") {
            tot <- apply(counts, 2, function(x) {
                stats::quantile(x[x > 0],
                    prob = UQ.pc
                )
            })
            if (sum(tot == 0) > 0) {
                stop(paste0(
                    "Cannot perform UQ normalization (percentile=",
                    UQ.pc, " ), not enough signal in sample(s) ",
                    paste(colnames(counts)[tot == 0], collapse = ", ")
                ))
            }
        } else {
            tot <- colSums(counts)
        }
        ncounts <- sweep(counts, 2, tot / stats::median(tot), "/")
    } else {
        ncounts <- counts
    }

    homolog.genes <- list()
    if (species != "hsapiens") {
        ncounts <- as.data.frame(ncounts)
        ncounts$human.gene.name <- rownames(ncounts)
        conversion.dict$human.gene.name <- rownames(conversion.dict)
        ncounts$id <- seq_len(nrow(ncounts))

        counts.transposed <- merge(ncounts, conversion.dict,
            by.x = "human.gene.name",
            all = FALSE, sort = FALSE
        )
        counts.transposed <- counts.transposed[order(counts.transposed$id), ]

        homolog.genes <- list(counts.transposed$Gene.name)

        counts.transposed$id <- NULL
        ncounts$id <- NULL
        ncounts$human.gene.name <- NULL
        ncounts <- data.matrix(ncounts)
        rm(counts.transposed)
    }
    nLR <- length(intersect(
        c(BulkSignalR:::.SignalR$BulkSignalR_LRdb$ligand,
        BulkSignalR:::.SignalR$BulkSignalR_LRdb$receptor),
        rownames(ncounts)
    ))
    if (nLR < min.LR.found) {
        stop(
            "Not enough LR genes (", nLR, " < ", min.LR.found,
            " were found).\n"
        )
    }

    bsrdm <- new("BSRDataModel",
        ncounts = ncounts, log.transformed = log.transformed,
        normalization = toupper(method), initial.organism = species,
        initial.orthologs = homolog.genes
    )
    new("SCSRNoNet",
        bsrdm.comp = as(bsrdm,"BSRDataModelComp"),
        populations = populations, paracrines = list(), 
        autocrines = list(), inf.param = list()
    )
} # SCSRNoNet


# Accessors & setters ========================================================

#' BSRDataModelComp object accessor
#'
#' @name bsrdmComp
#' @aliases bsrdmComp,SCSRNoNet-method
#' @param x SCSRNoNet object
#' @return bsrdmComp object
#' @export
#' @examples
#' print("bsrdmComp")
#' data(example_dataset, package = "SingleCellSignalR")
#' mat <- log1p(data.matrix(example_dataset[, -1])) / log(2)
#' rownames(mat) <- example_dataset[[1]]
#' rme <- rowMeans(mat)
#' mmat <- mat[rme > 0.05, ]
#' d <- dist(t(mmat))
#' h <- hclust(d, method = "ward.D")
#' pop <- paste0("pop_", cutree(h, 5))
#' scsrnn <- SCSRNoNet(mat,
#'     normalize = FALSE, method = "log-only",
#'     log.transformed = TRUE, populations = pop
#' )
#' bsrdmComp(scsrnn)
setMethod("bsrdmComp", "SCSRNoNet", function(x) x@bsrdm.comp)


#' BSRDataModelComp object setter (internal use only)
#'
#' @param x SCSRNoNet object
#' @param value value to be set for bsrdm.comp
#' @return returns \code{NULL}
#' @keywords internal
setMethod("bsrdmComp<-", "SCSRNoNet", function(x, value) {
    x@bsrdm.comp <- value
    methods::validObject(x)
    x
})


#' populations accessor
#'
#' @name populations
#' @aliases populations,SCSRNoNet-method
#' @param x SCSRNoNet object
#' @return populations
#' @export
setMethod("populations", "SCSRNoNet", function(x) x@populations)

#' populations setter (internal use only)
#'
#' @param x SCSRNoNet object
#' @param value value to be set for populations
#' @return returns \code{NULL}
#' @keywords internal
setMethod("populations<-", "SCSRNoNet", function(x, value) {
    x@populations <- value
    methods::validObject(x)
    x
})


#' paracrines accessor
#'
#' @name paracrines
#' @aliases paracrines,SCSRNoNet-method
#' @param x SCSRNoNet object
#' @return paracrines
#' @export
setMethod("paracrines", "SCSRNoNet", function(x) x@paracrines)


#' paracrines setter (internal use only)
#'
#' @param x SCSRNoNet object
#' @param value value to be set for paracrines
#' @return returns \code{NULL}
#' @keywords internal
setMethod("paracrines<-", "SCSRNoNet", function(x, value) {
    x@paracrines <- value
    methods::validObject(x)
    x
})


#' autocrines accessor
#'
#' @name autocrines
#' @aliases autocrines,SCSRNoNet-method
#' @param x SCSRNoNet object
#' @return autocrines
#' @export
setMethod("autocrines", "SCSRNoNet", function(x) x@autocrines)


#' autocrines setter (internal use only)
#'
#' @param x SCSRNoNet object
#' @param value value to be set for autocrines
#' @return returns \code{NULL}
#' @keywords internal
setMethod("autocrines<-", "SCSRNoNet", function(x, value) {
    x@autocrines <- value
    methods::validObject(x)
    x
})


# inference methods ============================================================

#' Inference of ligand-receptor interactions based on regulation
#'
#' This method computes all the autocrine and paracrine ligand-receptor
#' interactions. In addition, it is possible to restrict inferences to
#' autocrine or paracrine only, or do it for chosen cell populations.
#'
#' @name performInferences
#' @aliases performInferences,SCSRNoNet-method
#'
#' @param obj           A SCSRNoNet object.
#' @param autocrine     A logical indicating whether autocrine interactions
#'   should be inferred.
#' @param paracrine     A logical indicating whether paracrine interactions
#'   should be inferred.
#' @param selected.populations  A vector of cell population names to
#'   limit inferences to these very cell populations.
#' @param funDiffExpr   An optional function to compute the differential
#'   expression tables for each population. The function is called for
#'   one population at a time, with the \code{SCSRNoNet} object and
#'   the population name as parameters,
#'   and it must return a \code{\link[BulkSignalR]{BSRClusterComp-class}}
#'   object.
#' @param subsample.size  The number of cells to sample from a given
#'   population to estimate differential expression significance.
#'   If set to -1, then all the cells of each population are used thus
#'   implying that the same observation from pairs of populations harboring
#'   different sizes will result in potentially very different P-values.
#' @param n.resample  The number of times the sampling is performed.
#' @param check.pval.se  A logical indicating that an estimate of the standard
#'   error on the differential expression P-value must be computed, and
#'   the conservative estimate P-value + SE must replace the P-value estimate.
#' @param max.pval        The maximum P-value imposed to both the ligand
#'   and the receptor.
#' @param min.logFC       The minimum log2 fold-change allowed for
#'   both the receptor and the ligand.
#' @param neg.receptors     A logical indicating whether receptors are only
#'   allowed to be upregulated (FALSE), or up- and downregulated (TRUE).
#' @param min.LR.score   The minimum required LR-score.
#' @param fdr.proc      The procedure for adjusting P-values according to
#' \code{\link[multtest]{mt.rawp2adjp}}.
#' @param restrict.genes  A list of gene symbols that restricts ligands and
#'   receptors.
#' @param verbose  A logical activating reports on computation steps done.
#'
#' @details
#' The basis of the interaction inferences is the increased expression of
#' the ligand and the receptor in
#' their respective cell populations. To determine differential expression
#' it is necessary to compare cell populations and to generate data.frames
#' representing gene differential expression. The data.frame format
#' is defined in \code{\link[BulkSignalR]{BSRClusterComp-class}} class.
#' A default procedure is provided for generating all
#' the differential expression analyses, i.e.,
#' one per cell population comparing it to all the other cells. It relies
#' on Wilcoxon tests performed on a fixed number of cells per population
#' (\code{subsample.size} parameter) to avoid biases due to actual sizes. Such
#' sampling and significance analysis is performed \code{n.resample} times
#' and the median P-value is finally used.
#' It is possible to substitute a user-defined function for this purpose
#' in case a different notion of differential expression would be
#' preferred.
#'
#' Once all the individual cell population differential analyses are done,
#' autocrine and paracrine ligand-receptor interactions are generated.
#'
#' In addition to statistical significance estimated by taking the ligand
#' and the receptor P-values product, we compute SingleCellSignalR original
#' LR-score, based
#' on L and R cluster average expression.
#'
#' @return A SCSRNoNet with inferences set.
#'
#' @export
#'
#' @examples
#' data(example_dataset, package = "SingleCellSignalR")
#' print("performInferences")
#' mat <- log1p(data.matrix(example_dataset[, -1])) / log(2)
#' rownames(mat) <- example_dataset[[1]]
#' rme <- rowMeans(mat)
#' mmat <- mat[rme > 0.05, ]
#' d <- dist(t(mmat))
#' h <- hclust(d, method = "ward.D")
#' pop <- paste0("pop_", cutree(h, 5))
#' scsrnn <- SCSRNoNet(mat,
#'     normalize = FALSE, method = "log-only",
#'     log.transformed = TRUE, populations = pop
#' )
#' # performInferences
#' if (FALSE) {
#' scsrnn <- performInferences(scsrnn,
#'     verbose = TRUE, ,
#'     min.logFC = 1e-10, max.pval = 1)
#' }
#' @importFrom methods new
#' @importFrom BulkSignalR logTransformed mu 
#' @importFrom BulkSignalR ncounts BSRClusterComp addClusterComp
#' @importFrom BulkSignalR differentialStats comparison colClusterA colClusterB
#' @importFrom matrixStats rowMeans2
#' @importFrom foreach %do% %dopar%
#' @importFrom matrixTests row_wilcoxon_twosample
setMethod("performInferences", "SCSRNoNet", function(obj,
    autocrine = TRUE, paracrine = TRUE,
    selected.populations = NULL,
    funDiffExpr = NULL, subsample.size = 50,
    n.resample = 50, check.pval.se = FALSE,
    max.pval = 0.05, min.logFC = 1,
    min.LR.score = 0, neg.receptors = FALSE,
    restrict.genes = NULL,
    fdr.proc = c(
        "BH", "Bonferroni", "Holm", "Hochberg",
        "SidakSS", "SidakSD", "BY", "ABH", "TSBH"
    ),
    verbose = FALSE) {

    fdr.proc <- match.arg(fdr.proc)
    if (min.logFC <= 0) {
        stop("min.logFC must be >0")
    }
    if (!is.null(selected.populations) && 
        !all(selected.populations %in% populations(obj))) {
        stop("selected.populations contains unknown populations")
    }
    if (!is.null(funDiffExpr) && !is.function(funDiffExpr)) {
        stop("funDiffExpr must be a function")
    }
    if (!autocrine && !paracrine) {
        stop("One of paracrine or autocrine must be TRUE at least")
    }
    if (n.resample < 1){
        stop("n.resample must be > 0")
    }
    if ((subsample.size < 1) && (subsample.size != -1)){
        stop("subsample.size must be > 0 or equal to -1")
    }
    
    # store inference parameters
    inf.param <- list()
    inf.param$autocrine <- autocrine
    inf.param$paracrine <- paracrine
    if (is.null(selected.populations)) {
        selected.populations <- unique(populations(obj))
    }
    inf.param$selected.populations <- selected.populations
    inf.param$subsample.size <- subsample.size
    inf.param$n.resample <- n.resample
    mu <- mu(bsrdmComp(obj))

    # generate the differential tables and comparisons
    if (verbose) {
        message("Computing diffential expression tables:")
    }
    pop <- k <- NULL
    for (pop in selected.populations) {
        if (verbose) {
            message(" ", pop)
        }
        if (!is.null(funDiffExpr)) {
            bsrcc <- funDiffExpr(obj, pop)
            if (!is(bsrcc, "BSRClusterComp")) {
                stop("The provided funDiffExpr ",
                    "function has return a non BSRClusterComp object")
            }
        } else {
            A <- which(populations(obj) == pop)
            B <- which(populations(obj) != pop)
            if (subsample.size != -1){
                # resampling at subsample.size size
                d <- foreach::foreach(k = seq_len(n.resample),
                    .combine = cbind) %do% {
                    Ap <- sample(A, subsample.size, replace = TRUE)
                    Bp <- sample(B, subsample.size, replace = TRUE)
                    matrixTests::row_wilcoxon_twosample(
                        ncounts(bsrdmComp(obj))[, Ap],
                        ncounts(bsrdmComp(obj))[, Bp]
                    )$pvalue
                }
                diff <- apply(d, 1, stats::median, na.rm = TRUE)
                # NA typically caused ball all the values being equal
                diff[is.na(diff)] <- 1

                if (check.pval.se){
                    diff.se <- apply(d, 1, stats::mad, na.rm = TRUE)
                    diff.se[is.na(diff.se)] <- 1
                    diff <- diff + diff.se
                }
            }
            else{
                # no resampling, all the cells of each population are used
                diff <- matrixTests::row_wilcoxon_twosample(
                    ncounts(bsrdmComp(obj))[, A],
                    ncounts(bsrdmComp(obj))[, B]
                )$pvalue
            }
            
            expr <- matrixStats::rowMeans2(ncounts(bsrdmComp(obj))[, A])
            if (logTransformed(bsrdmComp(obj))) {
                logFC <- expr - matrixStats::rowMeans2(
                    ncounts(bsrdmComp(obj))[, B])
            } else {
                logFC <- (log1p(expr) - log1p(
                    matrixStats::rowMeans2(ncounts(
                        bsrdmComp(obj))[, B]))) / log(2)
            }
            tab <- data.frame(pval = diff, logFC = logFC, expr = expr)
            rownames(tab) <- rownames(ncounts(bsrdmComp(obj)))
            bsrcc <- BulkSignalR::BSRClusterComp(bsrdmComp(obj), A, B, tab)
        }
        bsrdmComp(obj) <- addClusterComp(
            bsrdmComp(obj),
            bsrcc, paste0(pop, "_vs_others")
        )
    }

    # autocrine inferences
    autocrines <- list()
    if (autocrine) {
        if (verbose) {
            message("Computing autocrine naive (network-free) interactions")
        }
        for (pop in selected.populations) {
            bsrcc <- BulkSignalR::comparison(bsrdmComp(obj))[[
                paste0(pop, "_vs_others")]]
            s <- differentialStats(bsrcc)
            naive.inter <- BulkSignalR:::.SignalR$BulkSignalR_LRdb[
            BulkSignalR:::.SignalR$BulkSignalR_LRdb$ligand %in% rownames(s) &
                BulkSignalR:::.SignalR$BulkSignalR_LRdb$receptor 
                %in% rownames(s), seq_len(2)]
            names(naive.inter) <- c("L", "R")
            if (!is.null(restrict.genes)) {
                naive.inter <- naive.inter[naive.inter$L %in% restrict.genes &
                    naive.inter$R %in% restrict.genes, ]
            }
            naive.inter$pval <- s[naive.inter$L, "pval"] * 
            s[naive.inter$R, "pval"]
            naive.inter$L.logFC <- s[naive.inter$L, "logFC"]
            naive.inter$R.logFC <- s[naive.inter$R, "logFC"]
            naive.inter$L.pval <- s[naive.inter$L, "pval"]
            naive.inter$R.pval <- s[naive.inter$R, "pval"]
            if (logTransformed(bsrdmComp(obj))) {
                sq <- s[naive.inter$L, "expr"] * s[naive.inter$R, "expr"]
            } else {
                sq <- log1p(s[naive.inter$L, "expr"]) / 
                log(2) * log1p(s[naive.inter$R, "expr"]) / log(2)
            }
            naive.inter$LR.score <- sq / (mu + sq)
            naive.inter <- naive.inter[naive.inter$L.pval <= max.pval &
                naive.inter$R.pval <= max.pval &
                naive.inter$L.logFC >= min.logFC &
                naive.inter$LR.score >= min.LR.score, ]
            if (neg.receptors) {
                naive.inter <- naive.inter[abs(
                    naive.inter$R.logFC) >= min.logFC, ]
            } else {
                naive.inter <- naive.inter[naive.inter$R.logFC >= min.logFC, ]
            }

            if (nrow(naive.inter) > 0) {
                autocrines <- c(autocrines, list(naive.inter))
                names(autocrines)[length(autocrines)] <- pop
            }
        }
    }

    # paracrine inferences
    paracrines <- list()
    if (paracrine) {
        if (verbose) {
            message("Computing paracrine naive (network-free) interactions")
        }
        for (source.pop in selected.populations) {
            src.bsrcc <- BulkSignalR::comparison(bsrdmComp(obj))[[paste0(
                source.pop, "_vs_others")]]
            src.s <- differentialStats(src.bsrcc)
            for (target.pop in selected.populations) {
                if (source.pop != target.pop) {
                    # if (verbose)
                    #   message("  from",source.pop,"to",target.pop)
                    bsrcc <- BulkSignalR::comparison(bsrdmComp(obj))[[paste0(
                        target.pop, "_vs_others")]]
                    tar.s <- differentialStats(bsrcc)
                    naive.inter <- BulkSignalR:::.SignalR$BulkSignalR_LRdb[
                    BulkSignalR:::.SignalR$BulkSignalR_LRdb$ligand 
                    %in% rownames(src.s) & 
                    BulkSignalR:::.SignalR$BulkSignalR_LRdb$receptor 
                    %in% rownames(tar.s),
                        seq_len(2)]
                    names(naive.inter) <- c("L", "R")
                    if (!is.null(restrict.genes)) {
                        naive.inter <- naive.inter[naive.inter$L %in%
                        restrict.genes & 
                        naive.inter$R %in% restrict.genes, ]
                    }
                    naive.inter$pval <- src.s[naive.inter$L, "pval"] * 
                    tar.s[naive.inter$R, "pval"]
                    naive.inter$L.logFC <- src.s[naive.inter$L, "logFC"]
                    naive.inter$R.logFC <- tar.s[naive.inter$R, "logFC"]
                    naive.inter$L.pval <- src.s[naive.inter$L, "pval"]
                    naive.inter$R.pval <- tar.s[naive.inter$R, "pval"]
                    if (logTransformed(bsrdmComp(obj))) {
                        sq <- sqrt(src.s[naive.inter$L, "expr"] 
                            * tar.s[naive.inter$R, "expr"])
                    } else {
                        sq <- sqrt(log1p(src.s[naive.inter$L, "expr"]) / 
                            log(2) * 
                            log1p(tar.s[naive.inter$R, "expr"]) /
                            log(2))
                    }
                    naive.inter$LR.score <- sq / (mu + sq)
                    naive.inter <- naive.inter[naive.inter$L.pval <= max.pval &
                        naive.inter$R.pval <= max.pval &
                        naive.inter$L.logFC >= min.logFC &
                        naive.inter$LR.score >= min.LR.score, ]
                    if (neg.receptors) {
                        naive.inter <- naive.inter[abs(naive.inter$R.logFC) 
                        >= min.logFC, ]
                    } else {
                        naive.inter <- naive.inter[naive.inter$R.logFC 
                        >= min.logFC, ]
                    }

                    if (nrow(naive.inter) > 0) {
                        paracrines <- c(paracrines, list(naive.inter))
                        names(paracrines)[
                        length(paracrines)] <- paste0(source.pop,
                            "_vs_", target.pop)
                    }
                }
            }
        }
    }

    # set the slots to contain the new inferences
    autocrines(obj) <- autocrines
    paracrines(obj) <- paracrines
    infParamSC(obj) <- inf.param

    obj
}) # performInferences


# convenience methods ==========================================================

#' Method to retrieve autocrine interactions
#'
#' @name getAutocrines
#' @aliases getAutocrines,SCSRNoNet-method
#'
#' @param obj           A SCSRNoNet object.
#' @param pop           Name of the cell population whose autocrine
#'   ligand-receptor interactions should be retrieved.
#'
#' @details
#' A  data.frame is returned that
#' contains the inferred interactions.
#'
#' @return A data.frame.
#'
#' @export
#'
#' @examples
#' print("getAutocrines")
#' data(example_dataset, package = "SingleCellSignalR")
#' mat <- log1p(data.matrix(example_dataset[, -1])) / log(2)
#' rownames(mat) <- example_dataset[[1]]
#' rme <- rowMeans(mat)
#' mmat <- mat[rme > 0.05, ]
#' d <- dist(t(mmat))
#' h <- hclust(d, method = "ward.D")
#' pop <- paste0("pop_", cutree(h, 5))
#' scsrnn <- SCSRNoNet(mat,
#'     normalize = FALSE, method = "log-only",
#'     log.transformed = TRUE, populations = pop
#' )
#' if(FALSE){
#' # infer ligand-receptor interactions from the comparison
#'
#' scsrnn <- performInferences(scsrnn,
#'     verbose = TRUE, ,
#'     min.logFC = 1e-10, max.pval = 1)
#' 
#' getAutocrines(scsrcn, "pop_1")
#' }
#' @importFrom methods is
setMethod("getAutocrines", "SCSRNoNet", function(obj, pop) {
    if (!is(obj, "SCSRNoNet")) {
        stop("obj must be of class SCSRNoNet")
    }
    if (!(pop %in% names(autocrines(obj)))) {
        stop("pop must be in the autocrine inferences")
    }

    autocrines(obj)[[pop]]
}) # getAutocrines


#' Method to retrieve paracrine interactions
#'
#' @name getParacrines
#' @aliases getParacrines,SCSRNoNet-method
#'
#' @param obj           A SCSRNoNet object.
#' @param source.pop    Name of the cell population from which
#'   ligand-receptor interactions originate (express the ligands).
#' @param target.pop    Name of the cell population towards which
#'   ligand-receptor interactions go (express the receptors).
#'
#' @details
#' A data.frame is returned that
#' contains the inferred interactions.
#'
#' @return A data.frame.
#'
#' @export
#'
#' @examples
#' print("getParacrines")
#' data(example_dataset, package = "SingleCellSignalR")
#' mat <- log1p(data.matrix(example_dataset[, -1])) / log(2)
#' rownames(mat) <- example_dataset[[1]]
#' rme <- rowMeans(mat)
#' mmat <- mat[rme > 0.05, ]
#' d <- dist(t(mmat))
#' h <- hclust(d, method = "ward.D")
#' pop <- paste0("pop_", cutree(h, 5))
#' scsrcn <- SCSRNoNet(mat,
#' normalize = FALSE, method = "log-only",
#' log.transformed = TRUE, populations = pop
#' )
#' if(FALSE){
#' # infer ligand-receptor interactions from the comparison
#'
#' scsrcn <- performInferences(scsrcn,
#'     verbose = TRUE, ,
#'     min.logFC = 1e-10, max.pval = 1)
#' 
#' getParacrines(scsrcn, "pop_1","pop_2")
#' }  
#' @importFrom methods is
setMethod("getParacrines", "SCSRNoNet", function(obj, source.pop,
    target.pop) {
    
    if (!is(obj, "SCSRNoNet")) {
        stop("obj must be of class SCSRNoNet")
    }
    pair <- paste0(source.pop, "_vs_", target.pop)
    if (!(pair %in% names(paracrines(obj)))) {
        stop(pair, " must be in the paracrine inferences")
    }

    paracrines(obj)[[pair]]
}) # getParacrines


#' Inference parameters accessor
#'
#' @name infParamSC
#' @aliases infParamSC,SCSRNoNet-method
#' @param x SCSRNoNet object.
#' @return infParamSC
#' @examples
#' data(example_dataset, package = "SingleCellSignalR")
#' mat <- log1p(data.matrix(example_dataset[, -1])) / log(2)
#' rownames(mat) <- example_dataset[[1]]
#' rme <- rowMeans(mat)
#' mmat <- mat[rme > 0.05, ]
#' d <- dist(t(mmat))
#' h <- hclust(d, method = "ward.D")
#' pop <- paste0("pop_", cutree(h, 5))
#' scsrcn <- SCSRNoNet(mat,
#' normalize = FALSE, method = "log-only",
#' log.transformed = TRUE, populations = pop
#' )
#' if(FALSE){
#' # infer ligand-receptor interactions from the comparison
#'
#' scsrcn <- performInferences(scsrcn,
#'     verbose = TRUE, ,
#'     min.logFC = 1e-10, max.pval = 1)
#' 
#' infParamSC(scsrcn)
#' } 
#' @export
setMethod("infParamSC", "SCSRNoNet", function(x) x@inf.param)


#' Inference parameters setter (internal use only)
#' @param x SCSRNoNet object.
#' @param value value to be set.
#' @return returns \code{NULL}
#' 
#' @keywords internal
setMethod("infParamSC<-", "SCSRNoNet", function(x, value) {
    x@inf.param <- value
    methods::validObject(x)
    x
})


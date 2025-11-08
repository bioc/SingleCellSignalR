#' SingleCellSignalR cellular network object
#'
#' An S4 class to represent data and inferences related to a single-cell
#' data set and aimed at predicting ligand-receptor interactions. Both
#' autocrine and paracrine interactions can be inferred and explored.
#'
#' @slot bsrdm.comp   A BulkSignalR BSRDataModelComp object containing
#' the expression matrix data as well as the comparisons between cell
#' populations.
#' @slot populations  A vector defining the cell population to which
#' each individual cell belongs to. The length of this vector must
#' be equal to the number of columns of the expression matrix.
#' @slot paracrines   A list of BulkSignalR BSRInferenceComp objects
#' containing all the paracrine ligand-receptor interactions between
#' all the possible pairs of cell populations.
#' @slot autocrines   A list of BulkSignalR BSRInferenceComp objects
#' containing all the autocrine ligand-receptor interactions for
#' each cell population.
#' @slot inf.param    Inference parameters.
#'
#' @details  This class is a container for all the data and inferences
#' related to a given single-cell project. Inferred interactions in
#' \code{paracrines} and \code{autocrines} can be further reduced to
#' eliminate redundancies such as multiple downstream pathways for the same
#' ligand-receptor interaction. See reduction functions
#' \code{"\link[=SCSRNet-class]{reduceToBestPathway}"},
#' \code{"\link[=SCSRNet-class]{reduceToLigand}"},
#' \code{"\link[=SCSRNet-class]{reduceToReceptor}"} and
#' \code{"\link[=SCSRNet-class]{reduceToPathway}"}.
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
#' scsrcn <- SCSRNet(mat,
#'     populations = pop, normalize = FALSE,
#'     method = "log-only", log.transformed = TRUE
#' )
#' @importFrom BulkSignalR coerce 
#' @importClassesFrom BulkSignalR BSRDataModel BSRDataModelComp
#' @importFrom methods is as
setClass("SCSRNet",
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
    "SCSRNet",
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
    "show", "SCSRNet",
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

#' Instantiate a SCSRNet object from expression data
#'
#' Take a matrix or a data frame containing single-cell RNA-sequencing
#' or proteomics as input data and return a SCSRNet
#' object ready for subsequent cellular network inference.
#' Normally, SCSRNet objects
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
#' @return A SCSRNet object with empty interactions.
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
#' @importFrom BulkSignalR coerce
#' @importFrom methods as
#' @examples
#' print("SCSRNet")
#' data(example_dataset, package = "SingleCellSignalR")
#' mat <- log1p(data.matrix(example_dataset[, -1])) / log(2)
#' rownames(mat) <- example_dataset[[1]]
#' rme <- rowMeans(mat)
#' mmat <- mat[rme > 0.05, ]
#' d <- dist(t(mmat))
#' h <- hclust(d, method = "ward.D")
#' pop <- paste0("pop_", cutree(h, 5))
#' scsrcn <- SCSRNet(mat,
#'     normalize = FALSE, method = "log-only",
#'     log.transformed = TRUE, populations = pop
#' )
#'
SCSRNet <- function(counts, populations, normalize = TRUE, 
    symbol.col = NULL, min.count = 1, prop = 0.1,
    method = c("UQ", "TC"), log.transformed = FALSE,
    min.LR.found = 80, species = "hsapiens", conversion.dict = NULL,
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
            stop("symbol.col must be the index ",
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
        stop("The read count matrix must be ",
            " provided with gene symbols as row names")
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
    new("SCSRNet",
        bsrdm.comp = as(bsrdm,"BSRDataModelComp"),
        populations = populations, paracrines = list(),
        autocrines = list(), inf.param = list()
    )
} # SCSRNet


# Accessors & setters ========================================================

setGeneric("bsrdmComp", signature="x",
    function(x) standardGeneric("bsrdmComp")
)
#' BSRDataModelComp object accessor
#'
#' @name bsrdmComp
#' @aliases bsrdmComp,SCSRNet-method
#' @param x SCSRNet object
#' @return bsrdmComp
#' @export
#' @examples
#' data(example_dataset, package = "SingleCellSignalR")
#' mat <- log1p(data.matrix(example_dataset[, -1])) / log(2)
#' rownames(mat) <- example_dataset[[1]]
#' rme <- rowMeans(mat)
#' mmat <- mat[rme > 0.05, ]
#' d <- dist(t(mmat))
#' h <- hclust(d, method = "ward.D")
#' clusters <- paste0("pop_", cutree(h, 5))
#'
#' scsrnn <- SCSRNet(mat,
#'     normalize = FALSE, method = "log-only", 
#'     min.count = 1, prop = 0.001,
#'     log.transformed = TRUE, populations = clusters
#' )
#'
#' bsrdmComp(scsrnn)
#'
setMethod("bsrdmComp", "SCSRNet", function(x) x@bsrdm.comp)


setGeneric("bsrdmComp<-", signature=c("x", "value"),
    function(x, value) standardGeneric("bsrdmComp<-")
)
#' BSRDataModelComp object setter (internal use only)
#' @param x SCSRNet object
#' @param value value to be set for bsrdm.comp
#' @return returns \code{NULL}
#' @keywords internal
setMethod("bsrdmComp<-", "SCSRNet", function(x, value) {
    x@bsrdm.comp <- value
    methods::validObject(x)
    x
})


setGeneric("populations", signature="x",
    function(x) standardGeneric("populations")
)
#' populations accessor
#'
#' @name populations
#' @aliases populations,SCSRNet-method
#' @param x SCSRNet object
#' @return populations
#' @examples
#' print("populations")
#' data(example_dataset, package = "SingleCellSignalR")
#' mat <- log1p(data.matrix(example_dataset[, -1])) / log(2)
#' rownames(mat) <- example_dataset[[1]]
#' rme <- rowMeans(mat)
#' mmat <- mat[rme > 0.05, ]
#' d <- dist(t(mmat))
#' h <- hclust(d, method = "ward.D")
#' pop <- paste0("pop_", cutree(h, 5))
#' scsrcn <- SCSRNet(mat,
#'     normalize = FALSE, method = "log-only",
#'     log.transformed = TRUE, populations = pop
#' )
#' 
#' populations(scsrcn)
#' @export
setMethod("populations", "SCSRNet", function(x) x@populations)


setGeneric("populations<-", signature=c("x", "value"),
    function(x, value) standardGeneric("populations<-")
)
#' populations setter (internal use only)
#' @param x SCSRNet object
#' @param value value to be set for populations
#' @return returns \code{NULL}
#' @keywords internal
setMethod("populations<-", "SCSRNet", function(x, value) {
    x@populations <- value
    methods::validObject(x)
    x
})


setGeneric("paracrines", signature="x",
    function(x) standardGeneric("paracrines")
)
#' paracrines accessor
#'
#' @name paracrines
#' @aliases paracrines,SCSRNet-method
#' @param x SCSRNet object
#' @return paracrines
#' @examples
#' print("paracrines")
#' data(example_dataset, package = "SingleCellSignalR")
#' mat <- log1p(data.matrix(example_dataset[, -1])) / log(2)
#' rownames(mat) <- example_dataset[[1]]
#' rme <- rowMeans(mat)
#' mmat <- mat[rme > 0.05, ]
#' d <- dist(t(mmat))
#' h <- hclust(d, method = "ward.D")
#' pop <- paste0("pop_", cutree(h, 5))
#' scsrcn <- SCSRNet(mat,
#'     normalize = FALSE, method = "log-only",
#'     log.transformed = TRUE, populations = pop
#' )
#' 
#' paracrines(scsrcn)
#' @export
setMethod("paracrines", "SCSRNet", function(x) x@paracrines)


setGeneric("paracrines<-", signature=c("x", "value"),
    function(x, value) standardGeneric("paracrines<-")
)
#' paracrines setter (internal use only)
#'
#' @param x SCSRNet object
#' @param value value to be set for paracrines
#' @return returns \code{NULL}
#' @keywords internal
setMethod("paracrines<-", "SCSRNet", function(x, value) {
    x@paracrines <- value
    methods::validObject(x)
    x
})


setGeneric("autocrines", signature="x",
    function(x) standardGeneric("autocrines")
)
#' autocrines accessor
#'
#' @name autocrines
#' @aliases autocrines,SCSRNet-method
#' @param x SCSRNet object
#' @return autocrines
#' @examples
#' print("autocrines")
#' data(example_dataset, package = "SingleCellSignalR")
#' mat <- log1p(data.matrix(example_dataset[, -1])) / log(2)
#' rownames(mat) <- example_dataset[[1]]
#' rme <- rowMeans(mat)
#' mmat <- mat[rme > 0.05, ]
#' d <- dist(t(mmat))
#' h <- hclust(d, method = "ward.D")
#' pop <- paste0("pop_", cutree(h, 5))
#' scsrcn <- SCSRNet(mat,
#'     normalize = FALSE, method = "log-only",
#'     log.transformed = TRUE, populations = pop
#' )
#' 
#' autocrines(scsrcn)
#' @export
setMethod("autocrines", "SCSRNet", function(x) x@autocrines)


setGeneric("autocrines<-", signature=c("x", "value"),
    function(x, value) standardGeneric("autocrines<-")
)
#' autocrines setter (internal use only)
#' @param x SCSRNet object
#' @param value value to be set for autocrines
#' @return returns \code{NULL}
#' @keywords internal
setMethod("autocrines<-", "SCSRNet", function(x, value) {
    x@autocrines <- value
    methods::validObject(x)
    x
})


setGeneric("infParamSC", signature="x",
    function(x) standardGeneric("infParamSC")
)
#' Inference parameters accessor
#'
#' @name infParamSC
#' @aliases infParamSC,SCSRNet-method
#' @param x SCSRNet object.
#' @export
#'
#' @examples
#' print("infParamSC")
#' data(example_dataset, package = "SingleCellSignalR")
#' mat <- log1p(data.matrix(example_dataset[, -1])) / log(2)
#' rownames(mat) <- example_dataset[[1]]
#' rme <- rowMeans(mat)
#' mmat <- mat[rme > 0.05, ]
#' d <- dist(t(mmat))
#' h <- hclust(d, method = "ward.D")
#' pop <- paste0("pop_", cutree(h, 5))
#' scsrcn <- SCSRNet(mat,
#'     normalize = FALSE, method = "log-only",
#'     log.transformed = TRUE, populations = pop
#' )
#' 
#' infParamSC(scsrcn)
#' @export
setMethod("infParamSC", "SCSRNet", function(x) x@inf.param)


setGeneric("infParamSC<-", signature=c("x", "value"),
    function(x, value) standardGeneric("infParamSC<-")
)
#' Inference parameters setter (internal use only)
#' @param x SCSRNet object.
#' @param value value to be set.
#' @return returns \code{NULL}
#' 
#' @keywords internal
setMethod("infParamSC<-", "SCSRNet", function(x, value) {
    x@inf.param <- value
    methods::validObject(x)
    x
})


# inference methods ============================================================

setGeneric("performInferences", signature="obj",
    function(obj, ...) standardGeneric("performInferences")
)
#' Inference of ligand-receptor interactions based on regulation
#'
#' This method computes all the autocrine and paracrine ligand-receptor
#' interactions. In addition, it is possible to restrict inferences to
#' autocrine or paracrine only, or do it for chosen cell populations.
#'
#' @name performInferences
#' @aliases performInferences,SCSRNet-method
#'
#' @param obj           A SCSRNet object.
#' @param autocrine     A logical indicating whether autocrine interactions
#'   should be inferred.
#' @param paracrine     A logical indicating whether paracrine interactions
#'   should be inferred.
#' @param selected.populations  A vector of cell population names to
#'   limit inferences to these very cell populations.
#' @param funDiffExpr   An optional function to compute the differential
#'   expression tables for each population. The function is called for
#'   one population at a time, with the \code{SCSRNet} object and
#'   the population name as parameters,
#'   and it must return a \code{\link[BulkSignalR]{BSRClusterComp-class}}
#'   object.
#' @param subsample.size  The number of cells to sample from a given
#'   population to estimate differential expression significance.
#' @param n.resample  The number of times the sampling is performed.
#' @param rank.p        A number between 0 and 1 defining the rank of the last
#'   considered target genes (see BulkSignalR documentation).
#' @param max.pval        The maximum P-value imposed to both the ligand
#'   and the receptor.
#' @param min.logFC       The minimum log2 fold-change allowed for
#'   both the receptor and the ligand.
#' @param neg.receptors     A logical indicating whether receptors are only
#'   allowed to be upregulated (FALSE), or up- and downregulated (TRUE).
#' @param fdr.proc      The procedure for adjusting P-values according to
#' \code{\link[multtest]{mt.rawp2adjp}}.
#' @param reference       Which pathway reference should be used ("REACTOME"
#'   for Reactome, "GOBP" for GO Biological Process,
#'   or "REACTOME-GOBP" for both).
#' @param max.pw.size     Maximum pathway size to consider from the pathway
#'   reference.
#' @param min.pw.size     Minimum pathway size to consider from the pathway
#'   reference.
#' @param min.positive    Minimum number of target genes to be found in a given
#'   pathway.
#' @param with.complex    A logical indicating whether receptor co-complex
#'   members should be included in the target genes.
#' @param min.t.logFC     The minimum log2 fold-change allowed for
#'   targets in case pos.targets or neg.targets are used.
#' @param pos.targets   A logical imposing that all the network targets must
#'   display positive logFC, i.e. logFC >= min.t.logFC.
#' @param neg.targets   A logical imposing that all the network targets must
#'   display negative logFC, i.e. logFC <= - min.t.logFC.
#' @param restrict.pw     A list of pathway IDs to restrict the application of
#'   the function.
#' @param restrict.genes  A list of gene symbols that restricts ligands and
#'   receptors.
#' @param use.full.network  A logical to avoid limiting the reference network
#' to the detected genes and use the whole reference network.
#' @param verbose  A logical activating reports on computation steps done.
#'
#' @details
#' The basis of the interaction inferences is the increased expression of
#' the ligand, the receptor, and the target genes below the receptor in
#' their respective cell populations. To determine differential expression
#' it is necessary to compare cell populations and to generate data.frames
#' representing gene differential expression. The data.frame format
#' is defined in \code{\link[BulkSignalR]{BSRClusterComp-class}} class.
#' A default procedure is provided for generating 
#' all the differential expression analyses, i.e.,
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
#' In addition to statistical significance estimated according to BulkSignalR
#' statistical model, we compute SingleCellSignalR original LR-score, based
#' on L and R cluster average expression.
#'
#' @return A SCSRNet object with inferences set.
#'
#' @export
#'
#' @examples
#' # prepare data
#' print("SCSRNet")
#' data(example_dataset, package = "SingleCellSignalR")
#' mat <- log1p(data.matrix(example_dataset[, -1])) / log(2)
#' rownames(mat) <- example_dataset[[1]]
#' rme <- rowMeans(mat)
#' mmat <- mat[rme > 0.05, ]
#' d <- dist(t(mmat))
#' h <- hclust(d, method = "ward.D")
#' pop <- paste0("pop_", cutree(h, 5))
#' scsrcn <- SCSRNet(mat,
#'     normalize = FALSE, method = "log-only",
#'     log.transformed = TRUE, populations = pop
#' )
#'
#' if (FALSE){
#' # infer ligand-receptor interactions from the comparison
#' scsrcn <- performInferences(scsrcn,
#'     verbose = TRUE, min.logFC = 1e-10,
#'     max.pval = 1
#' )
#' }
#' @importFrom methods new
#' @importFrom BulkSignalR logTransformed mu coerce 
#' @importFrom BulkSignalR ncounts BSRClusterComp addClusterComp
#' @importFrom BulkSignalR differentialStats comparison
#' @importFrom BulkSignalR colClusterA colClusterB updateInference
#' @importFrom foreach %do% %dopar%
#' @importFrom matrixStats rowMeans2
#' @importFrom matrixTests row_wilcoxon_twosample
setMethod("performInferences", "SCSRNet", function(obj, 
    autocrine = TRUE, paracrine = TRUE,
    selected.populations = NULL,
    funDiffExpr = NULL, subsample.size = 50,
    n.resample = 10, rank.p = 0.55,
    max.pval = 0.01, min.logFC = 1,
    min.LR.score = 0, neg.receptors = FALSE,
    pos.targets = FALSE, neg.targets = FALSE,
    min.t.logFC = 0.5, restrict.genes = NULL,
    use.full.network = FALSE,
    reference = c("REACTOME-GOBP", "REACTOME", "GOBP"),
    max.pw.size = 600, min.pw.size = 10, min.positive = 2,
    restrict.pw = NULL, with.complex = TRUE,
    fdr.proc = c(
        "BH", "Bonferroni", "Holm", "Hochberg",
        "SidakSS", "SidakSD", "BY", "ABH", "TSBH"),
    verbose = FALSE) {

    reference <- match.arg(reference)
    fdr.proc <- match.arg(fdr.proc)
    if (min.logFC <= 0) {
        stop("min.logFC must be >0")
    }
    if (min.t.logFC <= 0) {
        stop("min.t.logFC must be >0")
    }
    if (rank.p < 0 || rank.p > 1) {
        stop("rank.p must lie in [0;1]")
    }
    if (neg.targets && pos.targets) {
        stop("neg.targets and pos.targets cannot be TRUE simultaneously")
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

    # generate the differential tables and comparisons
    if (verbose) {
        message("Computing differential expression tables:")
    }
    pop <- k <- NULL
    for (pop in selected.populations) {
        if (verbose) {
            message(" ", pop)
        }
        if (!is.null(funDiffExpr)) {
            bsrcc <- funDiffExpr(obj, pop)
            if (!is(bsrcc, "BSRClusterComp")) {
                stop("The provided funDiffExpr function ",
                    "has return a non BSRClusterComp object")
            }
        } else {
            A <- which(populations(obj) == pop)
            B <- which(populations(obj) != pop)
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
            # typically caused by all the values being equals
            diff[is.na(diff)] <- 1 
            expr <- matrixStats::rowMeans2(ncounts(bsrdmComp(obj))[, A])
            if (logTransformed(bsrdmComp(obj))) {
                logFC <- expr - matrixStats::rowMeans2(
                    ncounts(bsrdmComp(obj))[, B])
            } else {
                logFC <- (log1p(expr) - log1p(matrixStats::rowMeans2(
                    ncounts(bsrdmComp(obj))[, B]))) / log(2)
            }
            tab <- data.frame(pval = diff, logFC = logFC, expr = expr)
            rownames(tab) <- rownames(ncounts(bsrdmComp(obj)))
            bsrcc <- BulkSignalR::BSRClusterComp(bsrdmComp(obj), A, B, tab)
        }
        bsrdmComp(obj) <- BulkSignalR::addClusterComp(
            bsrdmComp(obj),
            bsrcc, paste0(pop, "_vs_others")
        )
    }

    # universal inference
    pop <- selected.populations[1]
    universal.s <- BulkSignalR::differentialStats(
        BulkSignalR::comparison(bsrdmComp(obj))[[paste0(pop, "_vs_others")]])
    A <- BulkSignalR::colClusterA(
        BulkSignalR::comparison(bsrdmComp(obj))[[paste0(pop, "_vs_others")]])
    B <-  BulkSignalR::colClusterB(
        BulkSignalR::comparison(bsrdmComp(obj))[[paste0(pop, "_vs_others")]])
    universal.s$pval <- 1e-5
    universal.s$logFC <- 5
    universal.s$expr <- 5
    bsrcc <- BulkSignalR::BSRClusterComp(bsrdmComp(obj), A, B, universal.s)
    bsrdmComp(obj) <- BulkSignalR::addClusterComp(bsrdmComp(obj), 
    bsrcc, "scsr-universal")
    if (verbose) {
        message("Computing universal inference")
    }
    
    bsrinf.u <- BulkSignalR::BSRInferenceComp(bsrdmComp(obj),
        cmp.name = "scsr-universal",
        max.pval = 1, restrict.genes = restrict.genes,
        use.full.network = use.full.network,
        reference = reference, max.pw.size = max.pw.size,
        min.pw.size = min.pw.size, min.positive = min.positive,
        restrict.pw = restrict.pw, with.complex = with.complex
    )

    # autocrine inferences
    autocrines <- list()
    if (autocrine) {
        if (verbose) {
            message("Computing autocrine interactions")
        }
        for (pop in selected.populations) {
            bsrcc <- BulkSignalR::comparison(bsrdmComp(obj))[[
                paste0(pop, "_vs_others")]]
            bsrinf.comp <- BulkSignalR::updateInference(bsrinf.u,
                bsrcc, ncounts(bsrdmComp(obj)),
                rank.p = rank.p, max.pval = max.pval,
                min.logFC = min.logFC, min.LR.score = min.LR.score,
                neg.receptors = neg.receptors,
                pos.targets = pos.targets, neg.targets = neg.targets,
                min.t.logFC = min.t.logFC,
                min.positive = min.positive,
                fdr.proc = fdr.proc
            )
            if (is.null(bsrinf.comp)) {
                if (verbose) {
                    message("  No interaction selected for ", pop)
                }
            } else {
                autocrines <- c(autocrines, list(bsrinf.comp))
                names(autocrines)[length(autocrines)] <- pop
            }
        }
    }

    # paracrine inferences
    paracrines <- list()
    if (paracrine) {
        if (verbose) {
            message("Computing paracrine interactions")
        }
        for (source.pop in selected.populations) {
            src.bsrcc <- BulkSignalR::comparison(bsrdmComp(obj))[[
                paste0(source.pop,
                "_vs_others")]]
            for (target.pop in selected.populations) {
                if (source.pop != target.pop) {
                    bsrcc <- BulkSignalR::comparison(bsrdmComp(obj))[[
                        paste0(target.pop,
                        "_vs_others")]]
                    bsrinf.comp <- BulkSignalR::updateInference(
                        bsrinf.u, bsrcc,
                        ncounts(bsrdmComp(obj)), src.bsrcc = src.bsrcc,
                        rank.p = rank.p, max.pval = max.pval,
                        min.logFC = min.logFC, min.LR.score = min.LR.score,
                        neg.receptors = neg.receptors,
                        pos.targets = pos.targets, neg.targets = neg.targets,
                        min.t.logFC = min.t.logFC,
                        min.positive = min.positive,
                        fdr.proc = fdr.proc
                    )
                    if (is.null(bsrinf.comp)) {
                        if (verbose) {
                            message("  No interaction selected for ",
                                    source.pop, "_vs_", target.pop)
                        }
                    } else {
                        paracrines <- c(paracrines, list(bsrinf.comp))
                        names(paracrines)[length(paracrines)] <- paste0(
                          source.pop, "_vs_", target.pop)
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

setGeneric("getAutocrines", signature="obj",
    function(obj, ...) standardGeneric("getAutocrines")
)
#' Method to retrieve autocrine interactions
#'
#' @name getAutocrines
#' @aliases getAutocrines,SCSRNet-method
#'
#' @param obj           A SCSRNet object.
#' @param pop           Name of the cell population whose autocrine
#'   ligand-receptor interactions should be retrieved.
#'
#' @details
#' A  \code{\link[BulkSignalR]{BSRInferenceComp-class}} object is returned that
#' contains all the inferred (unfiltered) interactions.
#'
#' Interactions in tabular format can be obtained applying the
#' \code{\link[BulkSignalR]{LRinter}} accessor to the returned object. All
#' the BSRInferenceComp methods to reduce pathway, ligand, or receptor
#' redundancies can also be applied to the return object.
#' @return A BSRInferenceComp object.
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
#' scsrcn <- SCSRNet(mat,
#'     normalize = FALSE, method = "log-only",
#'     log.transformed = TRUE, populations = pop
#' )
#' if(FALSE){
#' # infer ligand-receptor interactions from the comparison
#'
#' scsrcn <- performInferences(scsrcn,
#' verbose = TRUE, min.logFC = 1e-10,
#' max.pval = 1)
#'
#' getAutocrines(scsrcn, "pop_1")
#' }
#' @importFrom methods is
setMethod("getAutocrines", "SCSRNet", function(obj, pop) {
    if (!is(obj, "SCSRNet")) {
        stop("obj must be of class SCSRNet")
    }
    if (!(pop %in% names(autocrines(obj)))) {
        stop("pop must be in the autocrine inferences")
    }

    autocrines(obj)[[pop]]
}) # getAutocrines


setGeneric("getParacrines", signature="obj",
    function(obj, ...) standardGeneric("getParacrines")
)
#' Method to retrieve paracrine interactions
#'
#' @name getParacrines
#' @aliases getParacrines,SCSRNet-method
#'
#' @param obj           A SCSRNet object.
#' @param source.pop    Name of the cell population from which
#'   ligand-receptor interactions originate (express the ligands).
#' @param target.pop    Name of the cell population towards which
#'   ligand-receptor interactions go (express the receptors and genes
#'   in pathways downstream the receptors).
#'
#' @details
#' A  \code{\link[BulkSignalR]{BSRInferenceComp-class}} object is returned that
#' contains all the inferred (unfiltered) interactions.
#'
#' Interactions in tabular format can be obtained applying the
#' \code{\link[BulkSignalR]{LRinter}} accessor to the returned object. All
#' the BSRInferenceComp methods to reduce pathway, ligand, or receptor
#' redundancies can also be applied to the return object.
#' @return A BSRInferenceComp object.
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
#' scsrcn <- SCSRNet(mat,
#' normalize = FALSE, method = "log-only",
#' log.transformed = TRUE, populations = pop
#' )
#' if(FALSE){
#' # infer ligand-receptor interactions from the comparison
#'
#' scsrcn <- performInferences(scsrcn,
#' verbose = TRUE, min.logFC = 1e-10,
#' max.pval = 1)
#'
#' getParacrines(scsrcn, "pop_1","pop_2")
#' }
#' @importFrom methods is
setMethod("getParacrines", "SCSRNet", function(obj, source.pop,
    target.pop) {
    if (!is(obj, "SCSRNet")) {
        stop("obj must be of class SCSRNet")
    }
    pair <- paste0(source.pop, "_vs_", target.pop)
    if (!(pair %in% names(paracrines(obj)))) {
        stop(pair, " must be in the paracrine inferences")
    }

    paracrines(obj)[[pair]]
}) # getParacrines

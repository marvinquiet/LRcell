#' Cell-type enrichment analysis for preranked gene set.
#'
#' This function wraps around \code{\link{LRcellCore}} in case of empty
#' inputs of the marker gene file and brain region.
#'
#' @name LRcell
#' @param gene.p Named vector of gene-level pvalues from DEG analysis, i.e.
#' DESeq2, LIMMA
#'
#' @param marker.g List of Cell-type specific marker genes derived from
#' single-cell RNA-seq. The name of the list is cell-type or cluster name, the
#' values are marker genes vectors or numeric named vectors. LRcell provides
#' marker genes list in different human/mouse brains, but users could use their
#' own marker gene list as input.
#' default: NULL
#'
#' @param species Either `mouse` or `human`, default: mouse.
#'
#' @param region Specific brain regions provided by LRcell. For mouse, LRcell
#' provides 9 brain regions: c("FC", "HC", "PC", "GP", "STR", "TH", "SN", "ENT",
#' "CB"). For human, LRcell provides c("pFC", "PBMC")
#'
#' @param method Either `logistic regression` or `linear regression`. Logistic
#' regression equally treats cell-type specific marker genes, however, if
#' certain values could determine the importance of marker genes, linear
#' regression can be performed, default: LR.
#'
#' @param min.size Minimal size of a marker gene set, will impact the balance of
#' labels
#'
#' @param sig.cutoff Cutoff for input genes pvalues, default: 0.05.
#'
#' @return A list with LRcell results. Each item represents a marker gene input.
#' Each item in this list is a statistics table. In the table, the row
#' represents the name of marker genes, and the columns are:
#'\itemize{
#'  \item ID The IDs of each marker genes, can be a cell type or cluster;
#'  \item genes_num How many marker genes are contributing to the analysis;
#'  \item coef The coefficients of Logistic Regression or Linear Regression;
#'  \item odds_ratio The odds ratio quantifies association in Logistic
#'  Regression;
#'  \item p-value The p-value calculated from the analysis;
#'  \item FDR The FDR after BH correction.
#'  \item lead_genes Genes that are contributing to the analysis;
#'}
#'
#' @import ExperimentHub
#' @import AnnotationHub
#'
#' @export
#' @examples
#' data(example_gene_pvals)
#' res <- LRcell(example_gene_pvals, species="mouse", region="FC", method="LR")
LRcell <- function(gene.p,
                marker.g = NULL,
                species = c("mouse", "human"),
                region = NULL,
                method = c("LR", "LiR"),
                min.size = 5,
                sig.cutoff = 0.05) {

    species <- match.arg(species)
    method <- match.arg(method)

    sc_marker_genes_list <- NULL
    regions <- NULL
    internal_data_ind <- TRUE # whether user uses internal data
    result_list <- NULL

    # user does not provide marker gene list
    if (is.null(marker.g)) {

        # user does not provide region information
        if (is.null(region)) {
            message("Because you did not choose a certain region, we will run all regions in this species for you.
It might take some time...")
            if (species == "mouse")
                regions <- MOUSE_REGIONS
            if (species == "human")
                regions <- HUMAN_REGIONS
        } else
            regions[1] <- region

        # get the list of marker genes
        for (reg in regions) {
            reg <- validate_region(species, reg)

            ## Use URL to acquire data, maybe needed in future if FTP server is ready
            # marker_genes_url <-
            #     "https://github.com/marvinquiet/LRcell/blob/master/marker_genes_lib/"
            # Github urls where enriched genes are stored
            # enriched_genes_url <- paste0(marker_genes_url,
            #                              paste0(species, '/', reg, 'enriched_genes.RDS'),
            #                              "?raw=true")
            #
            # N.TRIES <- 3  # try 3 times, refer to http://bioconductor.org/developers/how-to/web-query/
            # while (N.TRIES > 0L) {
            #     enriched_genes <- tryCatch(readRDS(url(enriched_genes_url)), error=identity)
            #     if (!inherits(enriched_genes, "error"))
            #         break
            #     N.TRIES <- N.TRIES - 1L
            # }
            #
            # if (N.TRIES == 0L) {
            #     stop("'getURL()' failed:",
            #          "\n  URL: ", enriched_genes_url,
            #          "\n  error: ", conditionMessage(enriched_genes),
            #          "\n Please check whether the Internect Connection is stable.")
            # }

            ## Use ExperimentHub to download data
            eh <- ExperimentHub::ExperimentHub()
            eh <- AnnotationHub::query(eh, "LRcellTypeMarkers")  ## query the data in AnnotationHub

            if (species == "mouse") {
                enriched_genes <- eh[[MOUSE_EXPHUB_MAPPING[reg]]]
            }

            if (species == "human") {
                enriched_genes <- eh[[HUMAN_EXPHUB_MAPPING[reg]]]
            }
            sc_marker_genes_list[[reg]] <- get_markergenes(enriched_genes, method)
        }
    } else {
        sc_marker_genes_list[["user"]] <- marker.g
        internal_data_ind <- FALSE
    }

    for (item in names(sc_marker_genes_list)) {
        marker_genes <- sc_marker_genes_list[[item]]
        res <- LRcellCore(gene.p = gene.p,
                        marker.g = marker_genes,
                        method = method,
                        min.size = min.size,
                        sig.cutoff = sig.cutoff)
        # if users use provided data, then add cell type to the dataframe, otherwise
        # LRcell uses ID as cell type
        if (internal_data_ind) {
            if (item %in% MOUSE_REGIONS) {
                ## set cell types for mouse brain data
                split_res <- strsplit(as.character(res$ID), "\\.")
                celltypes <- unlist(lapply(split_res, "[", 2))
                res$cell_type <- celltypes
            } else if ("pFC" == item) {
                ## set cell types for human pFC data
                split_res <- strsplit(as.character(res$ID), "_")
                celltypes <- unlist(lapply(split_res, "[", 1))
                res$cell_type <- celltypes
            } else if ("PBMC" == item) {
                ## set cell types for human PBMC data
                split_res <- strsplit(as.character(res$ID), "_")
                celltypes <- unlist(lapply(split_res, "[", 1))
                res$cell_type <- celltypes
            }
        } else {
            res$cell_type <- res$ID
        }

        result_list[[item]] <- res
    }
    result_list
}


#' Find most enriched cell types in bulk DE genes by Logistic Regression
#'
#' This is a function which takes marker genes from single-cell RNA-seq as
#' reference to calculate the enrichment of certain cell types in bulk DEG
#' analysis. We assume that bulk DEG is derived from certain cell-type specific
#' pattern.
#'
#' @name LRcellCore
#'
#' @param gene.p Named vector of gene-level pvalues from DEG analysis, i.e.
#' DESeq2, LIMMA
#'
#' @param marker.g List of Cell-type specific marker genes derived from
#' single-cell RNA-seq. The name of the list is cell-type or cluster name, the
#' values are marker genes vectors or numeric named vectors. LRcell provides
#' marker genes list in different human/mouse brains, but users could use their
#' own marker gene list as input.
#' default: NULL
#'
#' @param method Either `logistic regression` or `linear regression`. Logistic
#' regression equally treats cell-type specific marker genes, however, if
#' certain values could determine the importance of marker genes, linear
#' regression can be performed, default: LR.
#'
#' @param min.size Minimal size of a marker gene set, will impact the balance of
#' labels
#'
#' @param sig.cutoff Cutoff for input genes' pvalues, default: 0.05.
#'
#' @return A dataframe of LRcell statistics as described in \code{\link{LRcell}}.
#'
#' @import stats
#' @import utils
#'
#' @export
#'
#' @examples
#' data(mouse_FC_marker_genes)
#' data(example_gene_pvals)
#' res <- LRcellCore(example_gene_pvals, mouse_FC_marker_genes, method="LR")
LRcellCore <- function(gene.p,
                       marker.g,
                       method,
                       min.size = 5, sig.cutoff = 0.05) {

    if (typeof(marker.g) != "list")
        stop("Please make sure marker genes is a list with names indicating cell types or clusters.")

    if (method == "LiR") {
        # check whether marker genes are with statistics indicating significance
        if (all(lapply(marker.g, is_namedvector) == TRUE) == FALSE)
            stop("Please check each item in the list is a named vector with statistics for each gene.")
    }

    if (method == "LR") {
        # check whether marker genes are within each group
        if (all(lapply(marker.g, is_vector) == TRUE) == FALSE)
            stop("Please check each item is a vector containing gene symbols.")
    }

    # filter gene_pvals named vector
    gene.p <- gene.p[!is.na(gene.p)]
    gene.p[gene.p == 0] <- 10^(-15)
    gene_nlp <- -log(gene.p)

    # average duplicate genes
    if (length(names(gene_nlp)) != length(unique(names(gene_nlp)))) {
        dedup_gene_nlp <- tapply(gene_nlp, names(gene_nlp), mean)
        gene_nlp <- dedup_gene_nlp
    }
    sig_genes <- names(gene_nlp)[exp(-gene_nlp) < sig.cutoff]

    # filter single-cell marker gene list
    filtered_marker_genes <- marker.g[lapply(marker.g, length) >= min.size]

    # start analyze
    ids <- NA; coefs <- NA; pvals <- NA; genes_num <- NA; lead_genes <- NA
    idx <- 0;
    for (sub in names(filtered_marker_genes)) {
        marker_genes <- filtered_marker_genes[[sub]]
        if (method == "LiR") {
            # average duplicate marker genes
            dedup_marker_genes <- tapply(marker_genes, names(marker_genes), mean)
            marker_genes <- names(dedup_marker_genes)
        }
        matched_idx <- match(marker_genes, as.vector(names(gene_nlp), mode="character"))
        matched_idx <- matched_idx[!is.na(matched_idx)]
        matched_genes <- as.vector(names(gene_nlp), mode="character")[matched_idx]

        if (length(matched_genes) < min.size) next

        # start analyzing
        idx <- idx+1
        ids[idx] <- sub
        genes_num[idx] <- length(matched_genes)
        sig_DE_genes <- intersect(marker_genes, sig_genes)
        lead_genes[idx]<-paste(sig_DE_genes[order(sig_DE_genes)],collapse=", ")

        y <- rep(0, length(gene_nlp))
        if (method == "LR") {
            y[matched_idx] <- 1
            lr_df <- data.frame("y"=as.numeric(y), "x"=as.numeric(unname(gene_nlp)))
            # fit logistic regression
            glm.lr <- stats::glm(y ~ x,
                                   family = binomial(link="logit"),
                                   maxit = 100,
                                   lr_df)
            glm_res <- summary(glm.lr)
            coefs[idx] <- glm_res$coefficients["x","Estimate"]
            pvals[idx] <- glm_res$coefficients["x","Pr(>|z|)"]
        }

        if (method == "LiR") {
            matched_scores <- unname(dedup_marker_genes[matched_genes])
            y[matched_idx] <- as.numeric(matched_scores)
            lir_df <- data.frame("y"=y, "x"=as.numeric(unname(gene_nlp)))
            # fit linear regression
            glm.lir <- stats::glm(y ~ x,
                                family = gaussian(link="identity"),
                                maxit = 100,
                                lir_df)
            glm_res <- summary(glm.lir)
            coefs[idx] <- glm_res$coefficients["x","Estimate"]
            pvals[idx] <- glm_res$coefficients["x","Pr(>|t|)"]
        }
    }

    BH_fdr <- p.adjust(pvals,"BH")
    if (method == "LR") {
        # calculate the odds of marker gene being enriched when input pval
        # increase from 0.001 to 0.5
        odds_increase <- log(0.5)-log(0.001)
        odds_ratio <- exp(odds_increase*coefs)
        res <- data.frame("ID"=ids, "genes_num"=genes_num,
                        "coef"=coefs, "odds_ratio"=odds_ratio,
                        "p-value"=pvals, "FDR"=BH_fdr,
                        "lead_genes"=lead_genes)
    }

    if (method == "LiR") {
        res <- data.frame("ID"=ids, "genes_num"=genes_num,
                        "coef"=coefs, "p-value"=pvals, "FDR"=BH_fdr,
                        "lead_genes"=lead_genes)
    }
    res
}


#' Find most enriched cell types in bulk DE genes by Logistic Regression
#'
#' This is a function which takes marker genes from single-cell RNA-seq as
#' reference to calculate the enrichment of certain cell types in bulk DEG
#' analysis. This algorithm borrows from Marques et al, 2016 (https://science.sciencemag.org/content/352/6291/1326.long).
#'
#' @name LRcell_gene_enriched_scores
#'
#' @param expr Expression matrix with rows as genes and columns as cells, can be
#' an object of Matrix or dgCMatrix or a dataframe.
#'
#' @param annot Cell type annotation named vector with names as cell ids and
#' values as cell types.
#'
#' @param power The penalty on fraction of cells expressing the genes.
#'
#' @param parallel Whether to run it in parallel.
#'
#' @param n.cores How many cores to use in parallel mode.
#'
#' @return A numeric matrix with rows as genes and columns as cell types,
#' values are gene enrichment scores.
#'
#' @import BiocParallel
#' @export
LRcell_gene_enriched_scores <- function(expr,
                                     annot,
                                     power=1,
                                     parallel=TRUE,
                                     n.cores=4) {
    cat("Generate enrichment score for each gene..\n")
    gene_enriched_list <- NULL

    # check length
    if (ncol(expr) != length(annot))
        stop("Please check your provided cell type annotation is corresponding to the input expression dataframe.")
    # check corresponding cell names
    if(!all(sort(colnames(expr)) == sort(names(annot))))
        stop("Please check your provided cell type annotation is corresponding to the input expression dataframe.")


    if (parallel) { # do parallel
        BPPARAM <- BiocParallel::SnowParam(workers = n.cores,
                                           progressbar = TRUE,
                                           type = 'SOCK')

        gene_enriched_list <- BiocParallel::bplapply(rownames(expr),
                                                     enrich_posfrac_score,
                                                     expr,
                                                     annot,
                                                     power,
                                                     BPPARAM = BPPARAM)
    }
    else { # or just for-loop
        # set progress bar
        pb <- utils::txtProgressBar(max=nrow(expr), style=3)
        progress <- function(n) utils::setTxtProgressBar(pb, n)

        for (i in seq_len(nrow(expr))) {
            gene <- rownames(expr)[i]
            gene_enriched_list[[gene]] <- enrich_posfrac_score(gene,
                                                               expr,
                                                               annot,
                                                               power=power)
            progress(i)
        }
        close(pb)
    }

    cat("total gene number:", length(gene_enriched_list), '\n')
    enriched_genes <- do.call(rbind, gene_enriched_list)
    rownames(enriched_genes) <- rownames(expr)

    ## convert enriched gene scores dataframe into matrix
    numeric_enriched_genes <- apply(enriched_genes, c(1, 2), as.numeric)
    numeric_enriched_genes
}


#' Calculate enrichment scores for each cell type in a specific gene.
#'
#' This function takes a specific gene expression, cell type annotation and a
#' hyperparameter to calculate enrichment scores.
#'
#' @name enrich_posfrac_score
#'
#' @param gene Gene name from the expression matrix.
#'
#' @param expr Complete expression matrix with rows as genes and columns as cells.
#'
#' @param annot Cell type annotation named vector with names as cell ids and
#' values as cell types.
#'
#' @param power The penalty on fraction of cells expressing the genes
#'
#' @return Enrichment score list with cell type as names and enrichment score as
#' values.
enrich_posfrac_score <- function(gene, expr, annot, power=1) {
    gene.expr <- expr[gene, ]
    score_list <- list()
    for (celltype in unique(annot)) {
        cells <- names(annot[annot == celltype])
        cells_gene_expr <- gene.expr[cells]
        # enrichment in this cell type
        enrich <- (sum(cells_gene_expr)/length(cells_gene_expr)) / (sum(gene.expr)/length(gene.expr))
        # proportion of this gene express
        postfrac <- sum(cells_gene_expr > 0) / length(cells_gene_expr)
        score_list[[celltype]] <- enrich*(postfrac^power)
    }
    score_list
}


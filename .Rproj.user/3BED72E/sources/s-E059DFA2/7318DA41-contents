#' Wrapper function to run LRcell for preranked gene set.
#'
#' This function wraps around \link[LRcell]{LRcell_analysis} in case of empty
#' inputs of the marker gene file and brain region.
#'
#' @param gene_pvals Named vector of gene-level pvalues from DEG analysis, i.e.
#' DESeq2, LIMMA
#' @param sig_cutoff Cutoff for input genes' pvalues, default: 0.05.
#' @param sc_marker_genes List of Cell-type specific marker genes derived from
#' single-cell RNA-seq. The name of the list is cell-type or cluster name, the
#' values are marker genes vectors or numeric named vectors. LRcell provides
#' marker genes list in different human/mouse brains, but users could use their
#' own marker gene list as input.
#' default: NULL
#' @param species Either `mouse` or `human`, default: mouse.
#' @param region Specific brain regions provided by LRcell. For mouse, LRcell
#' provides 9 brain regions: c("FC", "HC", "PC", "GP", "STR", "TH", "SN", "ENT",
#' "CB"). For human, LRcell provides c("pFC")
#' @param method Either `logistic regression` or `linear regression`. Logistic
#' regression equally treats cell-type specific marker genes, however, if
#' certain values could determine the importance of marker genes, linear
#' regression can be performed, default: LR.
#' @param min_size Minimal size of a marker gene set, will impact the balance of
#' labels
#' @param plot Whether to generate the plot.
#' @return LRcell statistics.
#' @export
#' @examples
#' data(example_gene_pvals)
#' res <- LRcell(example_gene_pvals, region="FC")
LRcell <- function(gene_pvals,
                   sc_marker_genes = NULL,
                   sig_cutoff = 0.05,
                   species = c("mouse", "human"),
                   region = NULL,
                   method = c("LR", "LiR"),
                   min_size = 1,
                   plot=FALSE) {

    species <- match.arg(species)
    method <- match.arg(method)

    if (is.null(sc_marker_genes)) {
        if (is.null(region)) {
            message("Because you did not choose a certain region, we will run all regions in this specie for you. \n
                    It might take some time...")
            if (species == "mouse") {
                for (region in MOUSE_REGIONS) {
                    res <- LRcell_helper(gene_pvals, method, species, region, min_size, sig_cutoff)
                }
            } else if (species == "human") {
                for (region in HUMAN_REGIONS) {
                    res <- LRcell_helper(gene_pvals, method, species, region, min_size, sig_cutoff)
                }
            }
        } else {
            res <- LRcell_helper(gene_pvals, method, species, region, min_size, sig_cutoff)
        }
    }






}

#' Helper function to get single-cell marker genes
#'
#'  @param species `Mouse` or `Human`
#'  @param region Provided brain regions
#'  @return Provided single-cell marker genes list
LRcell_helper <- function(gene_pvals, method, species, region, min_size, sig_cutoff) {
    region <- validate_region(species, region)
    filepath <- system.file("extdata", paste0(species, '/', region, 'enriched_genes.RDS'),
                            package = "LRcell",
                            mustWork = TRUE)

    enriched_genes <- readRDS(filepath)
    sc_marker_genes <- get_markergenes(enriched_genes, method)

    res <- LRcell_analysis(gene_pvals,
                    sc_marker_genes, method,
                    min_size, sig_cutoff)
    res
}



#' Find most enriched cell types in bulk DE genes by Logistic Regression
#'
#' This is a function which takes marker genes from single-cell RNA-seq as
#' reference to calculate the enrichment of certain cell types in bulk DEG
#' analysis. We assume that bulk DEG is derived from certain cell-type specific
#' pattern.
#'
#' @export
#' @import stats
LRcell_analysis <- function(gene_pvals,
                   sc_marker_genes, method,
                   min_size = 1, sig_cutoff = 0.05) {

    if (typeof(sc_marker_genes) != "list")
        stop("Please make sure marker genes is a list with names indicating cell types or clusters.")

    if (method == "LiR") {
        # check whether marker genes are with statistics indicating significance
        if (all(lapply(sc_marker_genes, is_namedvector) == TRUE) == FALSE)
            stop("Please check each item in the list is a named vector with statistics for each gene.")
    }

    if (method == "LR") {
        # check whether marker genes are within each group
        if (all(lapply(sc_marker_genes, is_vector) == TRUE) == FALSE)
            stop("Please check each item is a vector containing gene symbols.")
    }

    # filter gene_pvals named vector
    gene_pvals <- gene_pvals[!is.na(gene_pvals)]
    gene_pvals[gene_pvals == 0] <- 10^(-15)
    gene_nlp <- -log(gene_pvals)

    # average duplicate genes
    if (length(names(gene_nlp)) != length(unique(names(gene_nlp)))) {
        dedup_gene_nlp <- tapply(gene_nlp, names(gene_nlp), mean)
        gene_nlp <- dedup_gene_nlp
    }
    sig_genes <- names(gene_nlp)[exp(-gene_nlp) < sig_cutoff]

    # filter single-cell marker gene list
    filtered_sc_marker_genes <- sc_marker_genes[sapply(sc_marker_genes, length) >= min_size]

    # start analyze
    ids <- NA; coefs <- NA; pvals <- NA; genes_num <- NA; lead_genes <- NA
    idx <- 0;
    for (sub in names(filtered_sc_marker_genes)) {
        marker_genes <- filtered_sc_marker_genes[[sub]]
        if (method == "LiR") {
            # average duplicate marker genes
            dedup_marker_genes <- tapply(marker_genes, names(maker_genes), mean)
            marker_genes <- names(dedup_marker_genes)
            marker_scores <- unname(dedup_marker_genes)
        }
        matched_idx <- match(marker_genes, as.vector(names(gene_pvals), mode="character"))
        matched_idx <- matched_idx[!is.na(matched_idx)]
        matched_genes <- as.vector(names(gene_pvals), mode="character")[matched_idx]

        if (length(matched_genes) < min_size) next

        # start analyzing
        idx <- idx+1
        ids[idx] <- sub
        genes_num[idx] <- length(matched_genes)
        sig_DE_genes <- intersect(marker_genes, sig_genes)
        lead_genes[idx]<-paste(sig_DE_genes[order(sig_DE_genes)],collapse=", ")

        y <- rep(0, length(gene_pvals))
        if (method == "LR") {
            y[matched_idx] <- 1
            lr_df <- data.frame("y"=as.numeric(y), "x"=as.numeric(unname(gene_nlp)))
            # fit logistic regression
            glm.lr <- stats::glm(y ~ x,
                                   family = binomial(link="logit"),
                                   lr_df)
            glm_res <- summary(glm.lr)
            coefs[idx] <- glm_res$coefficients["x","Estimate"]
            pvals[idx] <- glm_res$coefficients["x","Pr(>|z|)"]
        }

        if (method == "LiR") {
            matched_scores <- unname(marker_genes[matched_genes])
            y[matched_idx] <- as.numeric(matched_scores)
            lir_df <- data.frame("y"=y, "x"=as.numeric(unname(gene_nlp)))
            # fit linear regression
            glm.lir <- stats::glm(y ~ x,
                                  family = gaussian(link="identity"),
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



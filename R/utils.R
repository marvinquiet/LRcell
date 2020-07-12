#' mouse brain regions
MOUSE_REGIONS <- c("TH", "STR", "SN", "PC", "HC", "GP", "FC", "ENT", "CB")
HUMAN_REGIONS <- c("pFC")

#' Whether the input data is a named numeric vector
is_namedvector <- function(data) {
    is.vector(data) & is.numeric(data) & !is.null(names(data)) &
    !any(is.na(names(data)))
}

#' Whether the input data is a vector
is_vector <- function(data) {
    is.vector(data) & is.character(data) & is.null(names(data))
}

#' Whether the input region is valid
validate_region <- function(species, region) {
    if (species == "mouse")
        return(match.arg(region, MOUSE_REGIONS))
    else if (species == "human")
        return(match.arg(region, HUMAN_REGIONS))
}

#' Get top marker genes for each subcluster
get_markergenes <- function(enriched_genes, method=c("LR", "LiR"), topn=100) {
    method <- match.arg(method)

    marker_genes<- list()
    for (cluster in colnames(enriched_genes)) {
        enriched_values <- as.numeric(enriched_genes[, cluster])
        names(enriched_values) <- rownames(enriched_genes)
        sorted_enriched_genes <- sort(enriched_values, decreasing = TRUE)

        # For current version, we use an arbitrary number of genes for each cluster
        if (method == "LiR")
            marker_genes[[cluster]] <- head(sorted_enriched_genes, n=topn)
        else if (method == "LR")
            marker_genes[[cluster]] <- head(names(sorted_enriched_genes), n=topn)
    }
    marker_genes
}

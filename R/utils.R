#' mouse/human brain regions provided by LRcell package
#' @noRd
MOUSE_REGIONS <- c("TH", "STR", "SN", "PC", "HC", "GP", "FC", "ENT", "CB")
HUMAN_REGIONS <- c("pFC", "PBMC")

# Mapping data to ExperimentHub EH id
MOUSE_EXPHUB_MAPPING <- c("FC"="EH4548",
                    "CB"="EH4549",
                    "ENT"="EH4550",
                    "GP"="EH4551",
                    "PC"="EH4552",
                    "STR"="EH4553",
                    "SN"="EH4554",
                    "TH"="EH4555",
                    "HC"="EH4556")
HUMAN_EXPHUB_MAPPING <- c("pFC"="EH4557", "PBMC"="EH5420")


#' Whether the input data is a named numeric vector
#' @noRd
is_namedvector <- function(data) {
    is.vector(data) & is.numeric(data) & !is.null(names(data)) &
    !any(is.na(names(data)))
}

#' Whether the input data is a character vector
#' @noRd
is_vector <- function(data) {
    is.vector(data) & is.character(data) & is.null(names(data))
}

#' Whether the input region is valid
#' @noRd
validate_region <- function(species, region) {
    if (species == "mouse")
        return(match.arg(region, MOUSE_REGIONS))
    else if (species == "human")
        return(match.arg(region, HUMAN_REGIONS))
}

#' Get top marker genes for each subcluster
#'
#' @name get_markergenes
#'
#' @param enriched.g A return from \link[LRcell]{LRcell_gene_enriched_scores}
#' or from provided data
#'
#' @param method If LR, the return will be a list of genes; If LiR, the return
#' will be a list of named vector with names as genes and values as enriched
#' scores.
#'
#' @param topn Top N genes as marker genes.
#'
#' @return A list of top marker genes.
#'
#' @import utils
#' @export
#'
#' @examples
#' library(ExperimentHub)
#' eh <- ExperimentHub::ExperimentHub()
#' eh <- query(eh, "LRcellTypeMarkers")
#' # eh$title
#' enriched_genes <- eh[['EH4548']]
#' marker.g <- get_markergenes(enriched_genes, method="LR", topn=100)
get_markergenes <- function(enriched.g, method=c("LR", "LiR"), topn=100) {
    method <- match.arg(method)

    marker_genes<- list()
    for (cluster in colnames(enriched.g)) {
        enriched_values <- as.numeric(enriched.g[, cluster])
        names(enriched_values) <- rownames(enriched.g)
        sorted_enriched_genes <- sort(enriched_values, decreasing = TRUE)

        # For current version, we use an arbitrary number of genes for each cluster
        if (method == "LiR")
            marker_genes[[cluster]] <- head(sorted_enriched_genes, n=topn)
        else if (method == "LR")
            marker_genes[[cluster]] <- head(names(sorted_enriched_genes), n=topn)
    }
    marker_genes
}


# === only for curating cell types for provided data
#' Curate cluster name according to mouse cell types
#' @noRd
#' rename_clusters <- function(enriched_genes) {
#'     data(mouse_celltypes)
#'     clusters <- colnames(enriched_genes)
#'     clusters <- gsub('\\.', '-', clusters)
#'     celltypes_num <- table(mouse_celltypes[clusters])
#'     celltypes_cnt <- rep(1, length(celltypes_num))
#'     names(celltypes_cnt) <- names(celltypes_num)
#'
#'     new_clusternames <- rep(NA, length(clusters))
#'     for (i in seq_len(length(clusters))) {
#'         cluster <- clusters[i]
#'         celltype <- unname(mouse_celltypes[cluster])
#'         new_clusternames[i] <- paste(cluster,
#'                         paste(celltype, unname(celltypes_cnt[celltype]), sep='_'),
#'                         sep=".")
#'         celltypes_cnt[celltype] <- celltypes_cnt[celltype]+1
#'     }
#'
#'     new_enriched_genes <- enriched_genes
#'     colnames(new_enriched_genes) <- new_clusternames
#'     new_enriched_genes
#' }
#'
#' #' Change the cluster number into corresponding sub-cell types
#' #' @noRd
#' rename_by_extdata <- function() {
#'     files <- list.files(system.file("extdata", "mouse_bak", package = "LRcell"),
#'                         full.names = TRUE)
#'
#'     new_mouse_dir <- system.file("extdata", "new_mouse", package = "LRcell")
#'     for (enriched_file in files) {
#'         enriched_genes <- readRDS(enriched_file)
#'         new_enriched_genes <- rename_clusters(enriched_genes)
#'         saveRDS(new_enriched_genes, file=file.path(new_mouse_dir, basename(enriched_file)))
#'     }
#' }

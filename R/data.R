#' Mapping between subclusters and cell types in Mouse Brain
#'
#' A named vector containing the subclusters as name and cell types as values
#' in Mouse Brain. The cell types are pre-annotated by the dataset, which
#' includes: Endothelial, FibroblastLike, Mural, Oligodendrocytes,
#' Polydendrocytes, Astrocytes and Microglia.
#'
#' @format A named vector with 565 subclusters:
#' \describe{
#'     Named vector with name as subclusters and values as cell types.
#' }
#' @source \url{http://dropviz.org/} under tab `data`
"mouse_celltypes"

#' Example gene_pvals named vector from mouse experiment.
#'
#' A named vector containing gene symbols as name and p-values as values. This
#' is from a mouse Alzheimer's disease model (GEO: GSE90693), specifically 6
#' months after treatment in Frontal Cortex brain region. In this dataset, we
#' expect to see the Microglia as the most enriched cell type.
#'
#' @format A named vector with 23,420 items
#'
#' @source  \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90693}
#' `GSE90693_RawCountsData_TPR50_6months_AllRegions.txt.gz`
"example_gene_pvals"

#' Mouse/Human marker genes provided by package under `inst/extdata`
#'
#' A list containing gene uniqueness score for each cell type or cluster.
#'
#' @format A list {celtype: gene scores}
#'
#' @source \url{http://dropviz.org/} under tab `data`

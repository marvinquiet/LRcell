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

#' Example marker genes from mouse FC brain region.
#'
#' A list of marker genes with names indicating cell types. We selected top 100
#' enriched genes from each subcluster as marker genes list.
#'
#' @format A named vector with 81 subclusters in mouse Frontal Cortex:
#' \describe{
#'     Named vector with name as subclusters and values as marker genes.
#' }
#' @source \url{https://github.com/marvinquiet/LRcell/tree/master/marker_genes_lib}
"mouse_FC_marker_genes"


#' An example output of LRcell using data \link{example_gene_pvals} and
#' \link{mouse_FC_marker_genes}.
#'
#' @format A data frame with 81 rows as mouse FC sub-clusters and 8 variables:
#' \describe{
#'     \item{ID}{The IDs of each marker genes, can be a cell type or cluster}
#'     \item{genes_num}{How many marker genes are contributing to the analysis}
#'     \item{coef}{The coefficients of Logistic Regression or Linear Regression}
#'     \item{odds_ratio}{The odds ratio quantifies association in Logistic
#'     Regression}
#'     \item{p.value}{The p-value calculated from the analysis}
#'     \item{FDR}{The FDR after BH correction}
#'     \item{lead_genes}{Genes that are contributing to the analysis}
#'     \item{cell_type}{Cell typ ename}
#' }
#'
"example_LRcell_res"

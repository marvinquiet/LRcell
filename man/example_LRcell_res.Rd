% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{example_LRcell_res}
\alias{example_LRcell_res}
\title{An example output of LRcell using data \link{example_gene_pvals} and
\link{mouse_FC_marker_genes}.}
\format{
A data frame with 81 rows as mouse FC sub-clusters and 8 variables:
\describe{
    \item{ID}{The IDs of each marker genes, can be a cell type or cluster}
    \item{genes_num}{How many marker genes are contributing to the analysis}
    \item{coef}{The coefficients of Logistic Regression or Linear Regression}
    \item{odds_ratio}{The odds ratio quantifies association in Logistic
    Regression}
    \item{p.value}{The p-value calculated from the analysis}
    \item{FDR}{The FDR after BH correction}
    \item{lead_genes}{Genes that are contributing to the analysis}
    \item{cell_type}{Cell typ ename}
}
}
\usage{
data(example_LRcell_res)
}
\description{
An example output of LRcell using data \link{example_gene_pvals} and
\link{mouse_FC_marker_genes}.
}
\keyword{datasets}

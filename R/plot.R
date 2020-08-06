#' Manhattan plot for the enrichment of cell types
#'
#' This function draws out the LRcell result dataframe. In this function, we take
#' LRcell result dataframe and added cell types according to
#'
#' @name plot_manhattan_enrich
#' @usage plot_manhattan_enrich(lrcell_res)
#'
#' @param lrcell_res LRcell result dataframe
#'
#' @return A ggplot2 object
#'
#' @importFrom dplyr group_by arrange summarize
#' @import ggplot2
#' @import ggrepel
#'
#' @export
plot_manhattan_enrich <- function(lrcell_res,
                                  sig_cutoff = 0.05,
                                  label_topn = 5) {
  lr_tmp <- lrcell_res %>%
    dplyr::group_by(cell_type) %>%
    dplyr::arrange(cell_type)
  lr_tmp$pos <- as.numeric(rownames(lr_tmp))

  # axis center
  axisdf <- lr_tmp %>%
    dplyr::group_by(cell_type) %>%
    dplyr::summarize(center = (max(pos) + min(pos)) / 2)

  thres_topn <- ifelse(nrow(lr_tmp) > label_topn, label_topn, nrow(lr_tmp))
  repel_thres <- sort(lr_tmp$FDR)[thres_topn]
  g <- ggplot2::ggplot(lr_tmp, aes(x = pos, y = -log10(FDR))) +
    # add points
    ggplot2::geom_point(aes(color = as.factor(cell_type)), alpha = 0.8, size = 2) +
    ggplot2::geom_hline(yintercept = -log10(sig_cutoff), linetype = "dashed", color = "red") +
    # geom_label_repel(aes(label = ifelse(FDR <= repel_thres & FDR <= sig_cutoff, as.character(ID), "")),
    #   arrow = arrow(length = unit(0.03, "npc"), type = "closed", ends = "first"),
    #   nudge_y = 1
    # ) +
    ggrepel::geom_text_repel(aes(label = ifelse(FDR <= repel_thres & FDR <= sig_cutoff, as.character(ID), "")),
                    force = 5) +
    ggplot2::scale_x_continuous(label = axisdf$cell_type, breaks = axisdf$center) +
    # add labs
    ggplot2::labs(x = "clusters") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 1),
      panel.grid.minor.x = element_blank()
    )
  g
}

#' GSEA enrichment plot for the enrichment of different cell types
#'
#' This function draws out the LRcell result dataframe. In this function, we take
#' LRcell result dataframe and added cell types according to
#'
#' @name plot_gsea_enrich
#' @usage plot_gsea_enrich(gene_pvals, sc_marker_genes)
#'
#' @param gene.p Gene pvalues
#'
#' @param marker.g
#'
#' @return A ggplot2 object
#'
#' @import fgsea
#'
#' @export
plot_gsea_enrich <- function(gene.p, marker.g) {



}

#' Manhattan plot for the enrichment of cell types
#'
#' This function draws out the LRcell result dataframe. In this function, we take
#' LRcell result dataframe and added cell types according to
#'
#' @name plot_manhattan_enrich
#'
#' @param lrcell_res LRcell result dataframe.
#'
#' @param sig.cutoff The p-value cutoff showing significance result of LRcell.
#'
#' @param label.topn A numeric number showing how many significant cell types
#' will be labeled.
#'
#' @return A ggplot2 object
#'
#' @importFrom dplyr group_by arrange summarize
#' @import ggplot2
#' @import ggrepel
#' @import utils
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' data(example_LRcell_res)
#' plot_manhattan_enrich(example_LRcell_res)
plot_manhattan_enrich <- function(lrcell_res,
                                sig.cutoff = 0.05,
                                label.topn = 5) {
    lr_tmp <- lrcell_res %>%
        dplyr::group_by(.data$cell_type) %>%
        dplyr::arrange(.data$cell_type)
    lr_tmp$pos <- as.numeric(rownames(lr_tmp))

    # axis center
    axisdf <- lr_tmp %>%
        dplyr::group_by(.data$cell_type) %>%
        dplyr::summarize(center = (max(.data$pos) + min(.data$pos)) / 2)

    thres_topn <- ifelse(nrow(lr_tmp) > label.topn, label.topn, nrow(lr_tmp))
    repel_thres <- sort(lr_tmp$FDR)[thres_topn]
    g <- ggplot2::ggplot(lr_tmp, aes(x = .data$pos, y = -log10(.data$FDR))) +
        # add points
        ggplot2::geom_point(aes(color = as.factor(.data$cell_type)), alpha = 0.8, size = 2) +
        ggplot2::geom_hline(yintercept = -log10(sig.cutoff), linetype = "dashed", color = "red") +
        # geom_label_repel(aes(label = ifelse(FDR <= repel_thres & FDR <= sig.cutoff, as.character(ID), "")),
        #   arrow = arrow(length = unit(0.03, "npc"), type = "closed", ends = "first"),
        #   nudge_y = 1
        # ) +
        ggrepel::geom_text_repel(aes(label = ifelse(.data$FDR <= repel_thres & .data$FDR <= sig.cutoff, as.character(.data$ID), "")),
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

context("Test LRcell")

test_that("LRcell works on provided example dataset (Alzheimer's disease)", {
    data(example_gene_pvals)
    res <- LRcell(example_gene_pvals, region="FC")
    expect_equal(res$FC[which.min(res$FC$p.value), ]$cell_type, "Microglia")
})

test_that("LRcell enriched gene generation", {
    # generate a simulated gene*cell read counts matrix
    n.row <- 3; n.col <- 10
    sim.expr <- matrix(0, nrow=n.row, ncol=n.col)
    rownames(sim.expr) <- paste0("gene", 1:n.row)
    colnames(sim.expr) <- paste0("cell", 1:n.col)

    # generate a simulated annotation for cells
    sim.annot <- c(rep("celltype1", 3), rep("celltype2", 3), rep("celltype3", 4))
    names(sim.annot) <- colnames(sim.expr)

    sim.expr['gene1', ] <- c(3, 0, 2, 8, 10, 6, 1, 0, 0, 2) # marker gene for celltype2
    sim.expr['gene2', ] <- c(7, 5, 8, 1, 0, 5, 2, 3, 2, 1) # marker gene for celltype1
    sim.expr['gene3', ] <- c(8, 10, 6, 7, 8, 9, 5, 8, 6, 8) # house keeping

    # manual calculation on the result should be
    manual.enriched.df <- data.frame("celltype1" = c(0.3472222, 1.960784, 1.066667),
                              "celltype2" = c(2.5, 0.3921569, 1.066667),
                              "celltype3" = c(0.1171875, 0.5882353, 0.9))
    rownames(manual.enriched.df) <- rownames(sim.expr)

    # code result
    code.enriched.df <- LRcell_gene_enriched_scores(sim.expr, sim.annot, parallel = FALSE)

    expect_equal(manual.enriched.df['gene2', 'celltype2'],
                 unlist(code.enriched.df['gene2', 'celltype2']),
                 tolerance = .01)
    expect_equal(manual.enriched.df['gene1', 'celltype3'],
                 unlist(code.enriched.df['gene1', 'celltype3']),
                 tolerance = .01)
})

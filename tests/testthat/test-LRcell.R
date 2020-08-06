context("Test LRcell")

test_that("LRcell works on provided example dataset (Alzheimer's disease)", {
    data(example_gene_pvals)
    res <- LRcell(example_gene_pvals, region="FC")
    expect_equal(res$FC[which.min(res$FC$p.value), ]$cell_type, "Microglia")
})

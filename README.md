# LRcell tutorial

**Package version**: 0.99.0


## Introduction
The goal of LRcell is to identify specific sub-cell types that drives the changes 
observed in a bulk RNA-seq differential gene expression experiment. To achieve this, 
LRcell utilizes sets of cell marker genes acquired from single-cell RNA-sequencing (scRNA-seq) 
as indicators for various cell types in the tissue of interest. Next, for each cell type, 
using its marker genes as indicators, we apply Logistic Regression on the complete set of genes with 
differential expression p-values to calculate a cell-type significance p-value. 
Finally, these p-values are compared to predict which one(s) are likely to be responsible 
for the differential gene expression pattern observed in the bulk RNA-seq experiments. 
LRcell is inspired by the [LRpath](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2639007/) algorithm developed by Sartor et al., originally designed for pathway/gene set enrichment analysis. LRcell contains three major components: LRcell analysis, plot generation and marker gene selection. All modules in this package are written in R. 

**Pre-loaded marker genes**

LRcell provides marker genes in the Prefrontal Cortex (pFC) human brain region and 
nine mouse brain regions (Frontal Cortex, Cerebellum, Globus Pallidus, Hippocampus, Entopeduncular, 
Posterior Cortex, Striatum, Substantia Nigra and Thalamus). The human brain dataset describes 
single-cell gene expression profiles of control samples in [Major Depressive Disorder (MDD) disease studies](https://www.nature.com/articles/s41593-020-0621-y). The mouse brain dataset comes from 
a [whole mouse brain analysis](https://www.sciencedirect.com/science/article/pii/S0092867418309553) and the data is also available at [DropViz](http://dropviz.org/).

## Workflow

### Installation

The package can be installed by using `devtools`.
```{r}
install.packages("devtools")

library(devtools)
devtools::install_github("marvinquiet/LRcell")
```
Or you can download the code from the Github repo, build and install package yourself as follows:
```{bash}
git clone https://github.com/marvinquiet/LRcell.git
R CMD BUILD --no-build-vignettes LRcell  
R CMD INSTALL LRcell_1.0.0.tar.gz
```

To check whether LRcell package is successfully installed:
```{r}
library(LRcell)
```

### LRcell usage

Once we have LRcell package loaded, we can start using it to analyze the transcriptional
engagement of cell types or clusters. LRcell **takes both single-cell marker genes list 
and p-values of bulk DE genes as input** to calculate the enrichment of cell-type specific
marker genes in the ranked DE genes. 

As mentioned above, LRcell provides single-cell marker genes list in 1 human brain region
(Prefrontal Cortex) and 9 mouse brain regions (Frontal Cortex, Cerebellum, Globus Pallidus, 
Hippocampus, Entopeduncular, Posterior Cortex, Striatum, Substantia Nigra and Thalamus). 
The human data comes from control samples in [Major Depressive Disorder studies](https://www.nature.com/articles/s41593-020-0621-y). The mouse data comes from the [whole brain single-cell RNA-seq experiments](https://www.sciencedirect.com/science/article/pii/S0092867418309553). Another resource for this dataset is [DropViz](http://dropviz.org/).

Users can either download the enriched data from Github repo and process them into 
marker gene list or we can directly indicate which species and brain regions users 
want to use as the input.

#### Directly indicate species and brain region in LRcell

Compared to processing data yourself, a much easier way is to indicate species and brain regions. 
In this way, marker genes will be downloaded and loaded for you. Say, we want to use mouse Frontal Cortex
marker genes to do analysis on the example bulk experiment. (The example contains 23, 420 genes along with p-values calculated from [DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html). 
Data is processed from a mouse Alzheimer's disease model (GEO: [GSE90693](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90693)), 
which is 6 months after treatment in Frontal Cortex brain region.)

Here, we take Linear Regression as example:
```{r}
data("example_gene_pvals")
res <- LRcell(gene.p = example_gene_pvals,
              marker.g = NULL,
              species = "mouse",
              region = "FC",
              method = "LiR")
FC_res <- res$FC
# exclude leading genes for a better view
sub_FC_res <- subset(FC_res, select=-lead_genes)
head(sub_FC_res[order(sub_FC_res$p.value), ])
```

The top 6 significant lines of the result:

| ID                            | genes\_num | coef         | p\.value        | FDR             | cell\_type       |
|-------------------------------|------------|--------------|-----------------|-----------------|------------------|
| FC\_8\-2\.Astrocytes\_2       | 100        | 0\.023308388 | 3\.529518e\-124 | 2\.858909e\-122 | Astrocytes       |
| FC\_11\-4\.unknown\_3         | 100        | 0\.100232565 | 4\.384345e\-99  | 1\.775660e\-97  | unknown          |
| FC\_11\-1\.Microglia\_1       | 100        | 0\.055650878 | 1\.620703e\-82  | 4\.375899e\-81  | Microglia        |
| FC\_9\-1\.Oligodendrocytes\_1 | 100        | 0\.013475547 | 1\.478273e\-41  | 2\.993502e\-40  | Oligodendrocytes |
| FC\_9\-2\.Oligodendrocytes\_2 | 100        | 0\.032374837 | 1\.184612e\-30  | 1\.919071e\-29  | Oligodendrocytes |
| FC\_9\-3\.Oligodendrocytes\_3 | 100        | 0\.007531059 | 1\.066234e\-29  | 1\.439416e\-28  | Oligodendrocytes |


You can also try plotting out the result:
```{r} 
plot_manhattan_enrich(FC_res, sig.cutoff = .05, label.topn = 5)
```

#### Marker gene download and do LRcell analysis

Marker gene list downloading example (mouse, Frontal Cortex, Logistic Regression):
```{r}
github_url <- "https://github.com/marvinquiet/LRcell/blob/master/marker_genes_lib/"
enriched_genes_url <- paste0(github_url, "mouse/FCenriched_genes.RDS?raw=true")
enriched_genes <- readRDS(url(enriched_genes_url))
# get marker genes for LRcell in logistic regression
FC_marker_genes <- get_markergenes(enriched_genes, method="LR", topn=100)
```

Then, we can run LRcell analysis by using `LRcellCore()` function using Logistic Regression.
```{r}
res <- LRcellCore(gene.p = example_gene_pvals,
           marker.g = FC_marker_genes,
           method = "LR", min.size = 5, 
           sig.cutoff = 0.05, package.d = TRUE)
```
Here, `package.d` is an indicator showing whether users are using the LRcell-provided data. 
If you are using your own single-cell marker genes data, you an set this indicator as `package.d=FALSE`.


## Calculate gene enrichment scores from expression dataframe
LRcell also offers a function to calculate the enrichment scores of a read count expression matrix. 
The calculation is based on algorithm in this [science paper](https://science.sciencemag.org/content/352/6291/1326.long)
which takes enrichment of genes in certain cell types and fraction of cells expressing the gene into 
consideration. 

`LRcell_gene_enriched_scores()` function takes the gene-by-cell expression matrix and a cell-type annotation as input, which means a clustering algorithm and/or cell-type annotations should be done first. The columns of the expression matrix should match with cell names in the annotation vector. 

Below is a naive example on how to utilize this function to select marker genes.

**Example gene-by-cell expression matrix:**

| | cell1 | cell2 | cell3 | cell4 | cell5 | cell6 | cell7 | cell8| cell9 | cell10|
|---|---|---|---|---|---|---|---|---|---|---|
| gene1 | 3 | 0 | 2 | 8 | 10 | 6 | 1 | 0 | 0 | 2 |
| gene2 | 7 | 5 | 8 | 1 | 0 | 5 | 2 | 3 | 2 | 1 |
| gene3 | 8 | 10 | 6 | 7 | 8 | 9 | 5 | 8 | 6 | 8 |

**Example cell annotation**
| cell ID | cell1 |  cell2 | cell3 | cell4 | cell5 | cell6 | cell7 | cell8| cell9 | cell10|
|---|---|---|---|---|---|---|---|---|---|---|
| Annotation | celltype1 | celltype1 | celltype1 | celltype2 | celltype2 | celltype2 | celltype3 | celltype3 | celltype3 | celltype3 |

This toy example contains 3 genes and 10 cells. As you can tell from the matrix, __gene1__ is 
a marker gene of __celltype2__; __gene2__ is a marker gene of __celltype1__; __gene3__ is a house keeping gene.


```{r}
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

# generating the enrichment score 
enriched_res <- LRcell_gene_enriched_scores(expr = sim.expr,
                            annot = sim.annot, parallel = FALSE)
enriched_res
```

After using `LRcell_gene_enriched_scores()` function, we get the enriched scores for each gene in each cell type, which is:

| | celltype1 | celltype2 | celltype3|
|---|---|---|---|
| gene1 | 0.3472222 | 2.5 | 0.1171875 |
| gene2 | 1.960784 | 0.3921569 | 0.5882353 |
| gene3 | 1.066667 | 1.066667 | 0.9 |

For example, __gene2__ is a marker gene of __celltype1__, after calculation, it has the highest enrichment score in __celltype1__ among all three cell types.

Then we can choose marker genes using top 1 as threshold:

```{r}
marker_res <- get_markergenes(enriched.g = enriched_res,
                method = "LR", topn=1)
```

The result becomes:

| celltype | marker genes |
|---|---|
| celltype1 | gene2 |
| celltype2 | gene1 |
| celltype3 | gene3 |

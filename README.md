# LRcell tutorial

**Package version**: 1.0.0

Our paper is now on Briefings in Bioinformatics [here](https://doi.org/10.1093/bib/bbac063).


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

LRcell provides marker genes in the Prefrontal Cortex (pFC) human brain region, human PBMC and 
nine mouse brain regions (Frontal Cortex, Cerebellum, Globus Pallidus, Hippocampus, Entopeduncular, 
Posterior Cortex, Striatum, Substantia Nigra and Thalamus). The human brain dataset describes 
single-cell gene expression profiles of control samples in [Major Depressive Disorder (MDD) disease studies](https://www.nature.com/articles/s41593-020-0621-y). The mouse brain dataset comes from 
a [whole mouse brain analysis](https://www.sciencedirect.com/science/article/pii/S0092867418309553) and the data is also available at [DropViz](http://dropviz.org/).

## Workflow

### Installation

This is a R Bioconductor package and it can be installed by using `BiocManager`.
```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager") ## this will install the BiocManager package
BiocManager::install("LRcell")
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
(Prefrontal Cortex), 1 human PBMC and 9 mouse brain regions (Frontal Cortex, Cerebellum, Globus Pallidus, 
Hippocampus, Entopeduncular, Posterior Cortex, Striatum, Substantia Nigra and Thalamus). 

- The human brain data comes from control samples in [Major Depressive Disorder studies](https://www.nature.com/articles/s41593-020-0621-y). 
- The human PBMC data comes from volunteers in [HIV vaccine trial](https://www.biorxiv.org/content/10.1101/2020.10.12.335331v1) at time point of immediately before.
- The mouse data comes from the [whole brain single-cell RNA-seq experiments](https://www.sciencedirect.com/science/article/pii/S0092867418309553). Another resource for this dataset is [DropViz](http://dropviz.org/).

The data is stored in another Bioconductor ExperimentHub package named [LRcellTypeMarkers](https://github.com/marvinquiet/LRcellTypeMarkers). Users can access the data through ExperimentHub by:

```{r}
## for installing ExperimentHub
BiocManager::install("ExperimentHub")

## query data
library(ExperimentHub)
eh <- ExperimentHub::ExperimentHub()
eh <- AnnotationHub::query(eh, "LRcellTypeMarkers")
eh  ## this will list out EH number to access the calculated gene enrichment scores

## get mouse brain Frontal Cortex enriched genes
enriched.g <- eh[["EH4548"]]
marker.g <- get_markergenes(enriched_genes, method="LR", topn=100)
```

Users are also encouraged to process a read count matrix with cell annotation information into a gene enrichment scores matrix.
```{r eval=FALSE}
enriched.g <- LRcell_gene_enriched_scores(expr, annot, power=1, parallel=TRUE, n.cores=4)
```
Here, `expr` is a read count matrix with rows as genes and columns as cells. `annot` is a named-vector with names as cell names (which is in accordance with the column names of `expr`) and values as annotated cell types. `power` is a hyper-parameter controlling how much penalty for the proportion of cells expressing certain gene. `parallel` and `n.cores` are two parameters for executing function in parallel to accelerate the calculation.

#### Directly indicate species and brain region in LRcell

Compared to processing data yourself, a much easier way is to indicate species and brain region or tissue. 
In this way, marker genes are extracted from ExperimentHub accordingly. For example, we can use mouse Frontal Cortex
marker genes to do LRcell analysis on the example bulk experiment. (The example contains 23, 420 genes along with p-values calculated from [DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html). 
Data is processed from a mouse Alzheimer's disease model (GEO: [GSE90693](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90693)), 
which is 6 months after treatment in Frontal Cortex brain region.)

Here, we take method as Linear Regression as example:
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
| FC\_8\-2\.Astrocytes          | 100        | 0\.019701817 | 4\.761701e\-121 | 3\.856978e\-119 | Astrocytes       |
| FC\_11\-4\.unknown            | 90         | 0\.089774550 | 1\.327954e\-93  | 5\.378213e\-92  | unknown          |
| FC\_11\-1\.Microglia          | 98         | 0\.056507809 | 2\.718086e\-91  | 7\.338832e\-90  | Microglia        |
| FC\_9\-1\.Oligodendrocytes    | 100        | 0\.008969389 | 6\.032076e\-40  | 1\.221495e\-38  | Oligodendrocytes |
| FC\_9\-3\.Oligodendrocytes    | 98         | 0\.006969152 | 1\.660017e\-35  | 2\.689227e\-34  | Oligodendrocytes |
| FC\_9\-4\.Oligodendrocytes    | 97         | 0\.006447281 | 6\.246377e\-31  | 8\.432608e\-30  | Oligodendrocytes |


You can also try plotting out the result:
```{r} 
g <- plot_manhattan_enrich(FC_res, sig.cutoff = .05, label.topn = 5) ## a ggplot2 object is returned
g
```

## Calculate gene enrichment scores from expression dataframe
The gene enrichment scores calculation is based on algorithm in this [science paper](https://science.sciencemag.org/content/352/6291/1326.long)
which takes both enrichment of genes in certain cell types and fraction of cells expressing the gene into consideration. 

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

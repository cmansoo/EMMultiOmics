# Get Gene Groupings from Functional classification text file

## Description

Given a functional classification result of genes, the function divides genes into groups

## Usage

```r
get_grouping(gene_symbols, file_dir)
```

## Arguments

* `gene_symbols`: a vector of gene names as official symbols.
* `file_dir`: file path of functional classification obtained from DAVID tool. A txt file is expected. Please see vignette("gene_grouping") for guide

## Examples

```r
fp <- system.file("eg_fx_classification.txt", package="EMMultiOmics")
get_grouping(eg_gene_symbols, fp)
```


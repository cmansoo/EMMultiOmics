# Z-matrix builder

## Description

Builds Z matrix to be used with NEG_em

## Usage

```r
Zmat_builder(R2, G)
```

## Arguments

* `R2`: R-squared from EMVS result
* `G`: gene expression level

## Examples

```r
G <- GBM_data2$G
R2 <- GBM2_EMVS_res$R2

Zmat_builder(R2, G)
```


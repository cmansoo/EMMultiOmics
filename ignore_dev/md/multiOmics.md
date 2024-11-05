# multiOmics

## Description

Run first and second stage models

## Usage

```r
multiOmics(
  M,
  G,
  grouping,
  nu0 = 0.5,
  nu1 = 10^3,
  nu = 1,
  lambda = 1,
  a = 1,
  b = 1,
  EMVS_I = 10,
  EMVS_thresh = 1e-04,
  Y,
  C,
  Delta,
  a0,
  gstr,
  n_fold = 10,
  random_seed = NULL,
  NEG_I = 10,
  NEG_thresh = 1e-04
)
```

## Arguments

* `M`: DNA Methylation matrix
* `G`: gene expression level
* `grouping`: Gene grouping
* `nu0`: parameter 0 for spike-and-slab Gaussian mixture prior on $\beta$
* `nu1`: parameter 1 for spike-and-slab Gaussian mixture prior on $\beta$
* `lambda`: For the prior on $\sigma2$, an inverse gamma prior $\pi(\sigma2 | \gamma) = IG(\nu/2, \nu\lambda/2)$, <- what is lambda here?
* `a`: Parameter $a$ for the beta prior $\pi(\theta) \propto \theta a-1(1-\theta)b-1$, $a, b > 0$
* `b`: Parameter $b$ for the beta prior $\pi(\theta) \propto \theta a-1(1-\theta)b-1$, $a, b > 0$
* `EMVS_I`: Maximum number of iterations of EMVS
* `EMVS_thresh`: Threshold for convergence criterion for EMVS
* `Y`: clinical outcome
* `C`: clinical features
* `a0`: size sparsity
* `gstr`: prior, two options: "scale" or "1"
* `n_fold`: number of folds in K-fold cross validation
* `random_seed`: random seed for splitting folds for cross validation
* `NEG_I`: number of maximum iterations for NEG_em
* `NEG_thresh`: convergence criterion fpr NEG_em
* `Zmatrix`: Loading matrix
* `.mpmath`: function depends on mpmath package from python. pointer for mpmath package

## Details

`gstr` will take two options "scale" or "1." if `gstr` == "scale" then g = 1/N^2 where N = number of genes

## Examples

```r
# params
M <- GBM_data2$M
G <- GBM_data2$G
fp <- system.file("eg_fx_classification.txt", package="EMMultiOmics")
gene_grouping <- get_grouping(eg_gene_symbols, fp)
Y <- GBM_data2$Y
C <- GBM_data2$C
a0 <- 0.1
gstr <- "scale"
Delta <- GBM_data2$Delta
n_fold <- 10
random_seed <- 123
EMVS_I <- 3
NEG_I <- 3
EMVS_thresh <- 0.0001
NEG_thresh <- 0.0001

# run
multiOmics_obj <- multiOmics(
  M=M,
  G=G,
  grouping=gene_grouping,
  Y=Y,
  C=C,
  a0=a0,
  gstr=gstr,
  Delta=Delta,
  n_fold=n_fold,
  random_seed = random_seed,
  EMVS_I = EMVS_I,
  NEG_I=NEG_I,
  EMVS_thresh = EMVS_thresh,
  NEG_thresh = NEG_thresh
)

# plot
plot(multiOmics_obj)
```


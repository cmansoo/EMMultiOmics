# multiOmics_sensitivity

## Description

Run `multiOmics` for given sets of parameters `g` and `a`

## Usage

```r
multiOmics_sensitivity(
  gstr_vec,
  a0_vec,
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
  n_fold = 10,
  random_seed = NULL,
  NEG_I = 10,
  NEG_thresh = 1e-04
)
```

## Arguments

* `gstr_vec`: a vector of `g` values, two options: "scale" or 1
* `a0_vec`: a vector of `a` values
* `...`: Additional parameters for `multiOmics`

## Examples

```r
# params
M <- GBM_data2$M
G <- GBM_data2$G
fp <- system.file("eg_fx_classification.txt", package="EMMultiOmics")
gene_grouping <- get_grouping(eg_gene_symbols, fp)
Y <- GBM_data2$Y
C <- GBM_data2$C
Delta <- GBM_data2$Delta
n_fold <- 10
random_seed <- 123
EMVS_I <- 3
NEG_I <- 3
EMVS_thresh <- 0.0001
NEG_thresh <- 0.0001

a0 <- c(0.1, 1, 10, 50)
g <- list("scale", 1)


# run
multiOmics_s_obj <- multiOmics_sensitivity(
  M=M,
  G=G,
  grouping=gene_grouping,
  Y=Y,
  C=C,
  a0_vec=a0,
  gstr_vec=g,
  Delta=Delta,
  n_fold=n_fold,
  random_seed = random_seed,
  EMVS_I = EMVS_I,
  NEG_I=NEG_I,
  EMVS_thresh = EMVS_thresh,
  NEG_thresh = NEG_thresh
)

# summary
summary(multiOmics_s_obj)
plot(multiOmics_s_obj)
```


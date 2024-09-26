# 2nd stage modeling helper

## Description

Runs cross validation for 2nd stage modeling; outputs summary

## Usage

```r
cv_helper(Y, G, C, Delta, R2, a0, gstr, n_fold = 10, random_seed = NULL, ...)
```

## Arguments

* `Y`: clinical outcome
* `G`: gene expression level
* `C`: clinical features
* `Delta`: details...
* `R2`: vector of R2 (R-squared?) from first EMVS algorithm
* `a0`: details... (or this can be under `...`)
* `gstr`: details... (or this can be under `...`)
* `n_fold`: number of folds in K-fold cross validation
* `random_seed`: random seed for splitting folds for cross validation
* `...`: to supply any additional arguments to NEG_em, the arguments are `I` and `thresh`

## Details

... can take:

* `I` = number of iterations for `NEG_em`
* `thresh` = convergence criterion for `NEG_em`

## Examples

```r
G <- GBM_data2$G
Y <- GBM_data2$Y
C <- GBM_data2$C
Delta <- GBM_data2$Delta
R2 <- GBM2_EMVS_res$R2
a0 <- 0.1
gstr <- "scale"
n_fold <- 10
random_seed <- 123

cv_helper(Y, G, C, Delta, R2, a0, gstr, n_fold, random_seed,
          I=10, thresh=0.001)

# multiple a0
a0 <- c(0.1, 1, 10, 50)
lapply(a0, function(a) cv_helper(Y, G, C, Delta, R2, a0=a, gstr, n_fold, random_seed, 
                                 I=10, thresh=0.001)) |> 
  setNames(a0)
```


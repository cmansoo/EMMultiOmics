# Normal-Exponential-Gamma EM Algorithm

## Description

A short description...

## Usage

```r
NEG_em(
  Y,
  G,
  C,
  a0,
  gstr,
  Zmatrix,
  I = 10,
  thresh = 0.001,
  .mpmath = setup_mpmath()
)
```

## Arguments

* `Y`: clinical outcome
* `G`: gene expression level
* `C`: clinical features
* `a0`: details...
* `gstr`: details... (if "scale", then gstr = 1/N^2, where N=nrow(G))
* `Zmatrix`: details...
* `I`: number of maximum iterations
* `thresh`: convergence criterion
* `.mpmath`: function depends on mpmath package from python. pointer for mpmath package

## Examples

```r
G <- GBM_data2$G
Y <- GBM_data2$Y
C <- GBM_data2$C
R2 <- GBM2_EMVS_res$R2
Zmatrix <- Zmat_builder(R2, G)
a0 <- 0.1
gstr <- "scale"
mpmath <- setup_mpmath()

NEG_em(Y, G, C, a0, gstr, Zmatrix,
       I=10, thresh=0.001, .mpmath=mpmath)
```


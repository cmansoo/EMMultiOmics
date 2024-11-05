# Expectation Maximization Variable Selection Function

## Description

Performs EMVS algorithm for genomics data.
The default values for `nu0`, `nu1`, `lambda`, `a`, `b` were chosen per publication: *Veronika Ro<U+010D>kov<U+00E1> & Edward I. George (2013): EMVS: The EM Approach to Bayesian Variable Selection, Journal of the American Statistical Association*

## Usage

```r
EMVS(
  M,
  G,
  grouping,
  nu0 = 0.5,
  nu1 = 10^3,
  nu = 1,
  lambda = 1,
  a = 1,
  b = 1,
  I = 100,
  thresh = 1e-04
)
```

## Arguments

* `M`: DNA Methylation matrix
* `G`: Gene expression level
* `grouping`: Gene grouping
* `nu0`: parameter 0 for spike-and-slab Gaussian mixture prior on $\beta$
* `nu1`: parameter 1 for spike-and-slab Gaussian mixture prior on $\beta$
* `lambda`: For the prior on $\sigma2$, an inverse gamma prior $\pi(\sigma2 | \gamma) = IG(\nu/2, \nu\lambda/2)$, <- what is lambda here?
* `a`: Parameter $a$ for the beta prior $\pi(\theta) \propto \theta a-1(1-\theta)b-1$, $a, b > 0$
* `b`: Parameter $b$ for the beta prior $\pi(\theta) \propto \theta a-1(1-\theta)b-1$, $a, b > 0$
* `I`: Maximum number of iterations of EMVS
* `thresh`: Threshold for convergence criterion.

## Examples

```r
# example data
M <- GBM_data2$M
G <- GBM_data2$G
grouping <- GENE_GROUP2

EMVS(M, G, grouping, I=10, thresh=0.001)
# EMVS_par(M, G, grouping, I=10, thresh=0.001)
```


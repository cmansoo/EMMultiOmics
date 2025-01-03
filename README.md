This is the code and data used for the manuscript: Bayesian Shrinkage Models for Integration and Analysis of Multiplatform High-Dimensional Genomics Data, by Hao Xue, Sounak Chakraborty, and Tanujit Dey. Please contact Hao (hx222atcornell.edu) for any questions. 

Due to advance of biomedical technologies, it has become prevalent to collect biomedical data of the same patients from different platforms in clinical research, such as epigentic, gene expression and clinical features. As integrating genomic and clinical data from multiplatforms can provide complementary information for clinical study, it is indispensable to develop statistical tools to analyze multiomics and multiplatform data jointly. In this paper we propose a two-stage hierarchical Bayesian model which can integrate biomedical data from diverse platforms to find biomarkers associated with the clinical outcomes of interest. In addition to that, we have high dimensionality problem in several data platforms. In the first stage, we use Expectation Maximization based approach to learn the regulating mechanism between epigenetics (e.g. gene methylation) and gene expression while taking gene functional annotations into consideration. In the second stage, based on the regulating mechanism learned in the first stage, we employ modified Bayesian Lasso to select significant genes associated with clinical outcomes of interest while incorporating other available clinical features. Simulation studies suggest our model based data integration method shows higher accuracy in selecting predictive variables. Moreover, real data analysis based on a Giloblastoma (GBM) dataset reveals our method's potential to detect genes associated with GBM survival, as most of the selected genes are clinically and biomedically crucial in the pathology of GBM in existing scienti c literature.

This repo consists of:
(1) the main function: main.R, which is the main script used for real data analysis; 

(2) the ancillary functions: functions_final.R, which harbors all the functions used in (1), and lpcf.py, a function used in the second stage model；

(3) the data used in paper: GBM_data2.RData, which consists of gene expression and DNA methylation data of GBM patients collected by Network, C. G. A. R. et al. (2008). Comprehensive genomic characterization defines human glioblastoma genes and core pathways. Nature 455, 1061. and Brennan, C. W., Verhaak, R. G., McKenna, A., Campos, B., Noushmehr, H., Salama, S. R., Zheng, S., Chakravarty, D., Sanborn, J. Z., Berman, S. H., et al. (2013). The somatic genomic landscape of glioblastoma. Cell 155, 462-477.


## Installation
```
remotes::install_github("cmansoo/EMMultiOmics", build_vignettes=TRUE)
```

## Vignettes
```
library(EMMultiOmics)
vignette("EMMultiOmics")
vignette("gene_grouping")
```
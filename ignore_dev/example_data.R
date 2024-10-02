# install package
remotes::install_github("cmansoo/EMMultiOmics", build_vignettes=TRUE)

library(EMMultiOmics)
# 1. example data
GBM_data2
M <- GBM_data2$M # methylation
G <- GBM_data2$G # gene expression
C <- GBM_data2$C # clinical features
Y <- GBM_data2$Y # clinical outcome
R2 <- GBM2_EMVS_res$R2 # R2 result of EMVS
eg_gene_symbols # gene symbols used

# get gene grouping
fp <- system.file("eg_fx_classification.txt", package="EMMultiOmics")
shell.exec(fp) # example functional classification result from DAVID
GENE_GROUP <- get_grouping(eg_gene_symbols, fp) # gene grouping using DAVID


# 2. example running EMVS (runtime is very long with large I)
# the result of EMVS that I've run earlier is saved in an object `EMMultiOmics::GBM2_EMVS_res`
EMVS(M, G, grouping, I=10, thresh=0.001)

# 3. example running NEG_em
# setup
Zmatrix <- Zmat_builder(R2, G)
a0 <- 0.1 # a param
gstr <- "scale" # g param
mpmath <- setup_mpmath()

# run
NEG_result <- NEG_em(Y, G, C, a0, gstr, Zmatrix,
                     I=10, thresh=0.001, .mpmath=mpmath)

NEG_result
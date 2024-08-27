# setup
G <- GBM_data2$G
C <- GBM_data2$C
Y <- log(GBM_data2$Y)
M <- GBM_data2$M
DELTA <- GBM_data2$Delta
grouping <- GENE_GROUP2

# test 1
EMVS(M, G, grouping, I=2, THRESH=0.1)
# test 2
EMVS_par(M, G, grouping, I=2, THRESH=0.1)
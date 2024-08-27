my_test <- function(test, I=2, THRESH=0.1){
  # setup
  G <- GBM_data2$G
  C <- GBM_data2$C
  Y <- log(GBM_data2$Y)
  M <- GBM_data2$M
  DELTA <- GBM_data2$Delta
  grouping <- GENE_GROUP2
  
  if(test == 1) EMVS(M, G, grouping, I=I, THRESH=THRESH)
  else if(test == 2) EMVS_par(M, G, grouping, I=I, THRESH=THRESH)
}

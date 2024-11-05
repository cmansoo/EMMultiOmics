remotes::install_github("cmansoo/EMMultiOmics", build_vignettes=TRUE)
library(EMMultiOmics)

# setup for 1st stage
M <- GBM_data2$M
G <- GBM_data2$G
fp <- system.file("eg_fx_classification.txt", package="EMMultiOmics")
gene_grouping <- get_grouping(eg_gene_symbols, fp)
EMVS_result <- EMMultiOmics::GBM2_EMVS_res

# setup for 2nd stage
Y <- GBM_data2$Y
C <- GBM_data2$C
R2 <- EMVS_result$R2
Delta <- GBM_data2$Delta
.mpmath <- setup_mpmath()
a0 <- 0.1
gstr <- 1
n_fold <- 10
random_seed <- 123
N <- nrow(G)
K <- ncol(G)
J <- ncol(M)
L <- ncol(C)
# build Z matrix
Zmatrix <- Zmat_builder(R2, G)

# E-M loop step
# run cross validation
set.seed(random_seed)
folds <- .folds(1:N, num_fold=n_fold)

# initialize result objects
selected_biomarkers <- c()
res_cols <- c("n_selected_vars", "r2_train", "r2_test",
              "cindex_train", "cindex_test", "mse_train", "mse_test",
              colnames(C))

res_table <- data.frame(matrix(0, n_fold, length(res_cols))) |> 
  setNames(res_cols)

G_coeffs <- data.frame(matrix(0, n_fold, ncol(G))) |> 
  setNames(colnames(G))

# loop
for(jjj in 1:n_fold){
  test_indx <- folds[[jjj]]
  train_delta <- Delta[-test_indx]
  test_delta <- Delta[test_indx]
  
  NEG_output <- NEG_em(Y=Y[-test_indx],
                       G=as.matrix(G[-test_indx,]),
                       C=as.matrix(C[-test_indx,]),
                       a0=a0,
                       gstr=gstr,
                       Zmatrix=Zmatrix,
                       I = 50, # lower to make the computation faster
                       thresh=0.001,
                       .mpmath=.mpmath)
  
  # save result
  NEG_output$beta <- data.frame(NEG_output$beta)
  colnames(NEG_output$beta) <- c(colnames(G), colnames(C))
  estbeta <- NEG_output$beta[NEG_output$k, ]
  selected_biomarkers <- unique(c(selected_biomarkers,
                                  names(estbeta)[abs(estbeta) > 1e-5 & !names(estbeta) %in% colnames(C)]))
  # evaluate cross validation result
  pred_train <- cbind(G[-test_indx,], C[-test_indx,]) %*% as.numeric(estbeta)
  pred_test <- cbind(G[test_indx,], C[test_indx,]) %*% as.numeric(estbeta)
  res_table[jjj, "n_selected_vars"] <- sum(abs(estbeta) > 1e-5 & !names(estbeta) %in% colnames(C))
  res_table[jjj, "r2_train"] <- .rsq(pred_train, Y[-test_indx])
  res_table[jjj, "r2_test"] <- .rsq(pred_test, Y[test_indx])
  res_table[jjj, "cindex_train"] <- .cindx(pred_train, Y[-test_indx])
  res_table[jjj, "cindex_test"] <- .cindx(pred_test, Y[test_indx])
  res_table[jjj, "mse_train"] <- mean((pred_train - Y[-test_indx])^2)
  res_table[jjj, "mse_test"] <- mean((pred_test - Y[test_indx])^2)
  res_table[jjj, colnames(C)] <- estbeta[1, colnames(C)]
  G_coeffs[jjj, colnames(G)] <- estbeta[1, colnames(G)]
}

# take the average
performance_df <- data.frame(a=a0,
                             g=gstr,
                             as.list(colMeans(res_table)))
performance_df2 <- performance_df[, ! colnames(performance_df) %in% colnames(C)] |> 
  t() |> 
  `colnames<-`("value")

# assign gene groups from loading matrix
Z_df <- data.frame(names(gene_grouping), Zmatrix[,-1]) |> 
  `colnames<-`(c("gene_names", "M", "M_non_M", "non_M"))

selected_Z <- Z_df[Z_df$gene_names %in% selected_biomarkers, ]
M_effect <- selected_Z[selected_Z$M == 1, "gene_names"]
M_Mc_effect <- selected_Z[selected_Z$M_non_M == 1, "gene_names"]
Mc_effect <- selected_Z[selected_Z$non_M == 1, "gene_names"]
selected_genes <- list(`M Effect` = M_effect, `M + M^c Effect` = M_Mc_effect, `M^c Effect` = Mc_effect)


#### Calculate Highest posterior density interval
c_beta <- res_table[colnames(C)] |>
  colMeans() |> 
  as.matrix() |>
  `colnames<-`("beta")

g_beta <- G_coeffs[selected_biomarkers] |>
  colMeans() |>
  as.matrix()

# implement formula
CtC <- t(C) %*% C
vcov_mat <- solve(CtC + diag(x=1, nrow=dim(CtC)[1], ncol=dim(CtC)[2]))
std_err <- sqrt(diag(vcov_mat))
Y_hat <- Y - (G[, selected_biomarkers] %*% g_beta)
mu_beta_c <- vcov_mat %*% (t(C) %*% Y_hat) |> 
  `colnames<-`("mu_beta")

# sample from beta_c ~ MVN(mu_beta_c, vcov)
betas <- MASS::mvrnorm(n = 1e5, mu_beta_c, vcov_mat)
intv <- HDInterval::hdi(betas)
hpd_lower_95 <- intv[1,]
hpd_upper_95 <- intv[2,]

# output results
coef_result <- cbind(c_beta, mu_beta_c, std_err, hpd_lower_95, hpd_upper_95)
coef_result

# final_result <- list(
#   performance = performance_df2,
#   coeffs = coef_result,
#   selected_genes = selected_genes,
#   cv_result = res_table
# )
# 
# print(final_result)

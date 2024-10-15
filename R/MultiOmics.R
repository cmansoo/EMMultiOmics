#' multiOmics
#' 
#' Run first and second stage models
#' @param M DNA Methylation matrix
#' @param G Gene expression level
#' @param grouping Gene grouping
#' @param nu0 parameter 0 for spike-and-slab Gaussian mixture prior on \eqn{\beta}
#' @param nu1 parameter 1 for spike-and-slab Gaussian mixture prior on \eqn{\beta}
#' @param lambda For the prior on \eqn{\sigma2}, an inverse gamma prior \eqn{\pi(\sigma2 | \gamma) = IG(\nu/2, \nu\lambda/2)}, <- what is lambda here?
#' @param a Parameter \eqn{a} for the beta prior \eqn{\pi(\theta) \propto \theta a-1(1-\theta)b-1}, \eqn{a, b > 0}
#' @param b Parameter \eqn{b} for the beta prior \eqn{\pi(\theta) \propto \theta a-1(1-\theta)b-1}, \eqn{a, b > 0}
#' @param EMVS_I Maximum number of iterations of EMVS
#' @param EMVS_thresh Threshold for convergence criterion for EMVS
#' @param Y clinical outcome
#' @param G gene expression level
#' @param C clinical features
#' @param a0 size sparsity
#' @param gstr prior (if "scale", then gstr = 1/N^2, where N=nrow(G)) 
#' @param Zmatrix Loading matrix
#' @param NEG_I number of maximum iterations for NEG_em
#' @param NEG_thresh convergence criterion fpr NEG_em
#' @param .mpmath function depends on mpmath package from python. pointer for mpmath package
#' @param n_fold number of folds in K-fold cross validation
#' @param random_seed random seed for splitting folds for cross validation
#' 
#' @examples
#' # params
#' M <- GBM_data2$M
#' G <- GBM_data2$G
#' fp <- system.file("eg_fx_classification.txt", package="EMMultiOmics")
#' gene_grouping <- get_grouping(eg_gene_symbols, fp)
#' Y <- GBM_data2$Y
#' C <- GBM_data2$C
#' a0 <- 0.1
#' gstr <- "scale"
#' Delta <- GBM_data2$Delta
#' n_fold <- 10
#' random_seed <- 123
#' EMVS_I <- 3
#' NEG_I <- 3
#' EMVS_thresh <- 0.0001
#' NEG_thresh <- 0.0001
#' 
#' # run
#' multiOmics(
#'   M=M,
#'   G=G,
#'   grouping=grouping,
#'   Y=Y,
#'   C=C,
#'   a0=a0,
#'   gstr=gstr,
#'   Delta=Delta,
#'   n_fold=n_fold,
#'   random_seed = random_seed,
#'   EMVS_I = EMVS_I,
#'   NEG_I=NEG_I,
#'   EMVS_thresh = EMVS_thresh,
#'   NEG_thresh = NEG_thresh
#' )
#' @export
multiOmics <- function(
    M, G, grouping, nu0=0.5, nu1=10^3, nu=1, lambda=1, a=1, b=1, EMVS_I=10, EMVS_thresh=0.0001,
    Y, C, Delta, a0, gstr, n_fold=10, random_seed=NULL, NEG_I=10, NEG_thresh=0.0001
  ){
  
  # within MultiOmics, run EMVS
  EMVS_result <- EMVS(M,G, grouping=gene_grouping, I=EMVS_I, thresh=EMVS_thresh)

  # run 2nd stage
  N <- nrow(G)
  K <- ncol(G)
  J <- ncol(M)
  L <- ncol(C)
  .mpmath <- setup_mpmath()
  R2 <- EMVS_result$R2
  
  # build Z matrix
  Zmatrix <- Zmat_builder(R2, G)
  
  # E-M loop step
  # run cross validation
  if(is.null(random_seed)) set.seed(random_seed)
  folds <- caret::createFolds(1:N, k=n_fold)
  
  # initialize result objects
  selected_biomarkers <- c()
  res_cols <- c("n_selected_vars", "r2_train", "r2_test",
                "cindex_train", "cindex_test", "mse_train", "mse_test",
                colnames(C))
  
  res_table <- data.frame(matrix(0, n_fold, length(res_cols))) |> 
    setNames(res_cols)
  
  message("Starting 2nd stage.")
  message("Running cross validation...")
  pb <- txtProgressBar(min=0, max=n_fold, initial=0, style=3)
  
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
                         I=NEG_I,
                         thresh=NEG_thresh, 
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
    res_table[jjj, colnames(C)] <- estbeta[colnames(C)]
    setTxtProgressBar(pb, jjj)
  }
  
  # take the average
  performance_df <- data.frame(a=a0, g=ifelse(gstr==1/(N^2),'scale',gstr),
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
  
  # calculate SE?
  # reference https://stats.stackexchange.com/questions/44838/how-are-the-standard-errors-of-coefficients-calculated-in-a-regression/44841#44841
  c_beta <- estbeta[colnames(C)] |> as.matrix() |> t() |> `colnames<-`("beta")
  sigma_sq <- sum((Y - C %*% c_beta)^2) / (nrow(C) - ncol(C))
  vcov_mat <- sigma_sq * chol2inv(chol(t(C) %*% C)) 
  std_err <- sqrt(diag(vcov_mat))
  lower_95 <- c_beta - 1.96 * std_err
  colnames(lower_95) <- "lower_95"
  upper_95 <- c_beta + 1.96 * std_err
  colnames(upper_95) <- "upper_95"
  
  coef_result <- cbind(c_beta, std_err, lower_95, upper_95)
  
  final_result <- list(
    performance = performance_df2,
    coeffs = coef_result,
    selected_genes = selected_genes,
    cv_result = res_table
  )
  
  print(final_result[c("performance", "coeffs")])
  print(final_result$selected_genes)
  
  # return
  final_result
}



# 




#' summary graphs
#' 
#' Summary box plot of multiOmics model result
#' 
#' @param multiOmics_mod multiOmics model result
#' @export
summary_plot <- function(multiOmics_mod){
  res_table <- multiOmics_mod$cv_result
  unwanted_cols <- c("n_selected_vars", "r2_train", "r2_test",
                     "cindex_train", "cindex_test", "mse_train", "mse_test")
  res_table[, setdiff(names(res_table), unwanted_cols)] |>
    utils::stack() |> 
    setNames(c("beta", "clinical_x")) |> 
    ggplot2::ggplot(ggplot2::aes(x=beta)) +
    ggplot2::geom_boxplot() +
    ggplot2::facet_wrap(~clinical_x, scales = "free_x", ncol=1) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank())
  
}








# orphancode

#' 2nd stage modeling helper
#'
#' @description
#' Runs cross validation for 2nd stage modeling; outputs summary
#'
#'
#' @param Y clinical outcome
#' @param G gene expression level
#' @param C clinical features
#' @param Delta details...
#' @param R2 vector of R2 (R-squared?) from first EMVS algorithm
#' @param a0 details... (or this can be under `...`)
#' @param gstr details... (or this can be under `...`)
#' @param n_fold number of folds in K-fold cross validation
#' @param random_seed random seed for splitting folds for cross validation
#' @param ... to supply any additional arguments to NEG_em, the arguments are `I` and `thresh`
#'
#' @details
#' ... can take:
#' - `I` = number of iterations for `NEG_em`
#' - `thresh` = convergence criterion for `NEG_em`
#'
#' @examples
#' G <- GBM_data2$G
#' Y <- GBM_data2$Y
#' C <- GBM_data2$C
#' Delta <- GBM_data2$Delta
#' R2 <- GBM2_EMVS_res$R2
#' a0 <- 0.1
#' gstr <- "scale"
#' n_fold <- 10
#' random_seed <- 123
#'
#' cv_helper(Y, G, C, Delta, R2, a0, gstr, n_fold, random_seed,
#'           I=10, thresh=0.001)
#'
#' # multiple a0
#' a0 <- c(0.1, 1, 10, 50)
#' lapply(a0, function(a) cv_helper(Y, G, C, Delta, R2, a0=a, gstr, n_fold, random_seed,
#'                                  I=10, thresh=0.001)) |>
#'   setNames(a0)
#'
cv_helper <- function(Y, G, C, Delta, R2, a0, gstr, n_fold=10, random_seed=NULL, ...){
  ###### make this function into a summarizer next, but for now combine the summary and NEG em result,
  ###### NEG em result and summary function should really be separate in the end tho
  # check if the user has set up mpmath? or NEG_em is going to simply call this object by default?
  .mpmath <- setup_mpmath()
  
  N <- nrow(G) # user param should be G
  K <- ncol(G)
  Zmatrix <- Zmat_builder(R2, G)
  
  # E-M loop step
  if(is.null(random_seed)) set.seed(random_seed)
  
  folds <- caret::createFolds(1:N, k=n_fold)
  NEG_list <- list()
  
  # at the first loop Hao only used a0 = 0.1 but later he used c(0.1, 1, 10, 50) so this second stage function
  # should be loopable with diff params of a0 and gstr
  # a0 <- 0.1 # user param
  # gstr <- 1/N^2 # user param <- from the txt files Hao created gstr values he used were c(1, 4.162331e-05, "scale") which is the same as c(1, 1/N^2, 1/(N^2)) so scale i just the same as 1/N^2
  # but not sure where to find these from the code he wrote main.R
  
  message("Starting cross validation...")
  pb <- txtProgressBar(min=0, max=n_fold, initial=0, style=3)
  
  # for(jjj in 1:n_fold){
  #   test_indx <- folds[[jjj]]
  #   train_delta <- Delta[-test_indx]
  #   test_delta <- Delta[test_indx]
  #   lst <- NEG_em(Y=Y[-test_indx],
  #                 G=as.matrix(G[-test_indx,]),
  #                 C=as.matrix(C[-test_indx,]),
  #                 a0=a0,
  #                 gstr=gstr,
  #                 Zmatrix=Zmatrix,
  #                 # I=10,
  #                 # thresh=0.001,
  #                 ...,
  #                 .mpmath=.mpmath)
  #   # save result
  #   NEG_list[[jjj]] <- lst
  #
  #   # track progress
  #   setTxtProgressBar(pb, jjj)
  #   if(jjj == 1) order_suffix <- "st"
  #   else if(jjj == 2) order_suffix <- "nd"
  #   else if(jjj == 3) order_suffix <- "rd"
  #   else order_suffix <- "th"
  #   msg_str <- paste0("\n", jjj, order_suffix, " fold complete")
  #
  #   message(msg_str)
  # }
  #
  # NEG list saved now
  # return
  # NEG_list
  
  # now create the summary
  # initialize result table
  cols <- c("a", "g", "size_support",
            "r2_train", "r2_test", "cindex_train", "cindex_test",
            "mse_train", "mse_test", "AGE", "PRIOR_GLIOMA",
            "SEX", "PRETREATMENT_HISTORY") # does this need to be parameterized? not sure what hes doing here
  final_table <- data.frame(matrix(0, nrow=8, ncol=length(cols))) # where is the nrow and ncol from?
  colnames(final_table) <- cols
  # final_table
  box_tables <- data.frame()
  # box_tables
  
  ## some looping?
  # result table
  ppp <- 1
  for(.gstr in gstr){
    for (.a0 in a0){
      res_cols <- c("size_support", "r2_train", "r2_test", "cindex_train", "cindex_test",
                    "mse_train", "mse_test",
                    # covariates that should be specific to the dataset? - ask Hao
                    "AGE", "PRIOR_GLIOMA", "SEX", "PRETREATMENT_HISTORY")
      res_table <- data.frame(matrix(0, n_fold, length(res_cols))) |>
        setNames(res_cols)
      
      selected_biomarkers <- c()
      
      for(jjj in 1:n_fold){
        
        test_indx <- folds[[jjj]]
        train_delta <- Delta[-test_indx]
        test_delta <- Delta[test_indx]
        lst <- NEG_em(Y=Y[-test_indx],
                      G=as.matrix(G[-test_indx,]),
                      C=as.matrix(C[-test_indx,]),
                      a0=a0,
                      gstr=gstr,
                      Zmatrix=Zmatrix,
                      # I=10,
                      # thresh=0.001,
                      ...,
                      .mpmath=.mpmath)
        # save result
        NEG_list[[jjj]] <- lst
        
        # track progress
        setTxtProgressBar(pb, jjj)
        if(jjj == 1) order_suffix <- "st"
        else if(jjj == 2) order_suffix <- "nd"
        else if(jjj == 3) order_suffix <- "rd"
        else order_suffix <- "th"
        msg_str <- paste0("\n", jjj, order_suffix, " fold complete")
        
        message(msg_str)
        
        
        # test_indx = folds[[jjj]]
        output = NEG_list[[jjj]]
        output$beta = data.frame(output$beta)
        colnames(output$beta) = colnames(G)
        estbeta = output$beta[output$k,]
        selected_biomarkers = unique(c(selected_biomarkers, names(estbeta)[abs(estbeta)> 1e-5]))
        pred_train = cbind(G[-test_indx,], C[-test_indx,]) %*% as.numeric(estbeta)
        pred_test = cbind(G[test_indx,], C[test_indx,]) %*% as.numeric(estbeta)
        res_table[jjj, 'size support'] = sum(abs(estbeta) > 1e-5)
        res_table[jjj,'r2_train'] = .rsq(pred_train, Y[-test_indx] )
        res_table[jjj,'r2_test']= .rsq(pred_test,Y[test_indx])
        res_table[jjj,'cindex_train']= .cindx(pred_train,Y[-test_indx])
        res_table[jjj,'cindex_test']= .cindx(pred_test,Y[test_indx])
        res_table[jjj,'mse_train']= mean((pred_train-Y[-test_indx])^2)
        res_table[jjj,'mse_test']= mean((pred_test-Y[test_indx])^2)
        res_table[jjj,c("AGE","PRIOR_GLIOMA","SEX","PRETREATMENT_HISTORY")] = estbeta[1, 1001:1004] # 1001:1004 from where?
      }
      final_table[ppp,1]=.a0
      final_table[ppp,2]=ifelse(.gstr==1/(N^2),'scale',.gstr)
      final_table[ppp,3:length(cols)]=round(colMeans(res_table),3) # 3:13 from where?
      box_table=data.frame(a=.a0,g=ifelse(.gstr==1/(N^2),'scale',.gstr),
                           AGE=res_table$AGE,
                           PRIOR_GLIOMA=res_table$PRIOR_GLIOMA,
                           SEX=res_table$SEX,
                           PRETREATMENT_HISTORY=res_table$PRETREATMENT_HISTORY)
      box_tables = rbind(box_tables,box_table)
      ppp=ppp+1
    }
  }
  box_tables$group = paste0('g=',box_tables$g,',a=',box_tables$a)
  # return
  list(final_table=final_table, box_table=box_tables)
}

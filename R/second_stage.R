#' Normal-Exponential-Gamma EM Algorithm
#' 
#' @description
#' A short description...
#' 
#' 
#' @param Y clinical outcome
#' @param G gene expression level
#' @param C clinical features
#' @param a0 details...
#' @param gstr details... (if "scale", then gstr = 1/N^2, where N=nrow(G)) 
#' @param Zmatrix details...
#' @param I number of maximum iterations
#' @param thresh convergence criterion
#' @param .mpmath function depends on mpmath package from python. pointer for mpmath package
#' 
#' @examples
#' G <- GBM_data2$G
#' Y <- GBM_data2$Y
#' C <- GBM_data2$C
#' R2 <- GBM2_EMVS_res$R2
#' Zmatrix <- Zmat_builder(R2, G)
#' a0 <- 0.1
#' gstr <- "scale"
#' mpmath <- setup_mpmath()
#' 
#' NEG_em(Y, G, C, a0, gstr, Zmatrix,
#'        I=10, thresh=0.001, .mpmath=mpmath)
#' @export
NEG_em <- function(Y, G, C, a0, gstr, Zmatrix, I=10, thresh=0.001, .mpmath=setup_mpmath()){ # dont need mypath, and simply output? can we have an option to ouput in terminal?
  # setup mpmath here? or is there a way to set it up upon installation of the package?
  # can use `source(...)` for this purpose? e.g. within the script define the function and 
  # source the function call inside a function?
  # if cannot automate, add a step in NEG_em, or lpcf that checks for mpmath pacakge
  # tell the user to install python and mpmath as it is a dependency
  # tell them to run .setup_mpmath() function as well (if the user has to run, then no longer private)
  
  # setup params
  N <- nrow(G)
  K0 <- ncol(G)
  L <- ncol(C)
  K <- K0 + L
  # I <- 300
  a <- a0
  g <- if(gstr == "scale") 1/N^2 else as.numeric(gstr)
  
  alpha <- 1 # these below and alpha should be parameterized or no?
  gam <- 1 #
  s1 <- 1 #
  s2 <- 2 #
  .c <- 1 #
  d <- 1 #
  
  # Store the values of expectations in E-step
  Elambdaj <- as.vector(rep(0,K0))
  Elambdaj2 <- as.vector(rep(0,K0)) 
  nz <- length(Zmatrix[1,])
  # initial values
  sigmak <- rep(1,I)
  betak <- matrix(0,nrow=I,ncol=K)
  bbk <- matrix(1,nrow=I,ncol=nz)
  v1 <- as.vector(rep((-(2*a+1+s1)),K0))
  v2 <- as.vector(rep((-(2*a+1+s2)),K0))
  vb <- as.vector(rep(-2*(a+0.5),K0))
  
  # starting iteration
  for(kkk in 2:I){
    #E-step
    
    #### should this be message instead?
    # print(bbk[kkk-1,])
    
    
    hf=1/((Zmatrix)%*%bbk[kkk-1,])
    ZZ=as.vector(abs(betak[kkk-1,1:K0])*sqrt(hf)/(sqrt(g)*sigmak[kkk-1]))
    lPCFb <- lpcf(K0,vb,ZZ, .mpmath)
    lPCF1 <- lpcf(K0,v1,ZZ, .mpmath)
    lPCF2 <- lpcf(K0,v2,ZZ, .mpmath)
    logmar=as.vector(rep(0,K0))
    for (h in 1:K0){
      logmar[h]=as.numeric(log(a*(2^a)*sqrt(hf[h])/(sqrt(pi)*sqrt(g)*sigmak[kkk-1]))+lgamma(a+0.5)+(betak[kkk-1,h])^2*hf[h]/(4*(g)*sigmak[kkk-1]^2)+lPCFb[h])
    }
    for(h in 1:K0){
      Elambdaj[h]=as.numeric(exp(lgamma(2*a+s1+1)+(0.5*s1+0.5)*log(hf[h])-log(sqrt(g)*sigmak[kkk-1]*gamma(a)*2^(a+0.5*s1))+betak[kkk-1,h]^2*(hf[h])/(4*(g)*sigmak[kkk-1]^2)-logmar[h]+lPCF1[h]))
      Elambdaj2[h]=as.numeric(exp(lgamma(2*a+s2+1)+(0.5*s2+0.5)*log(hf[h])-log(sqrt(g)*sigmak[kkk-1]*gamma(a)*2^(a+0.5*s2))+betak[kkk-1,h]^2*(hf[h])/(4*(g)*sigmak[kkk-1]^2)-logmar[h]+lPCF2[h]))
    }
    
    #M-step
    #update Beta through Adaptive LASSO, Zou (2006) (using "lqa" package)
    Bayelsso_gen = glmnet::glmnet(x=as.matrix(G),
                          y=Y-as.vector(as.matrix(C)%*%betak[kkk-1,(K0+1):K]),
                          intercept=FALSE, standardize=FALSE, alpha=1,
                          family="gaussian", lambda=2*sqrt(2)/sqrt(g)*sigmak[kkk-1]*(sum(Elambdaj[1:K0])/N)/(2*N),
                          penalty.factor =as.vector(Elambdaj[1:K0]))
    betak[kkk,1:K0]=coef(Bayelsso_gen)[-1]

    Bayelsso_clinical = glmnet::glmnet(x=as.matrix(C),
                               y=Y-as.vector(as.matrix(G)%*%betak[kkk,1:K0]),
                               intercept=FALSE, standardize=FALSE, alpha=0,
                               family="gaussian", lambda=1/N)
    betak[kkk,(K0+1):K]=coef(Bayelsso_clinical)[-1]
    
    ### plot for what?
    plot(betak[kkk,],main=kkk)
    
    
    #update sigma (correct the typo on denominator in the paper)
    sigmak[kkk]=as.numeric(((sqrt(2)*sum(abs(Elambdaj[1:K0]/sqrt(g)*betak[kkk,1:K0]))+sqrt(2*(sum(abs(Elambdaj[1:K0]/sqrt(g)*betak[kkk,1:K0])))^2+4*((t(Y-cbind(G,C)%*%betak[kkk,])%*%(Y-cbind(G,C)%*%betak[kkk,]))^1+betak[kkk,(K0:K)]%*%betak[kkk,(K0:K)]+2*d)*(N+K+2*.c+2+L)))/(2*(N+K+2*.c+2+L))))
    
    #update b
    Q2=function(bk){
      sum(-a*log(1/(t(bk)%*%t(Zmatrix)))-Elambdaj2[1:K0]*(t(bk)%*%t(Zmatrix)))+sum((alpha-1)*log(bk)-gam*bk)
    }
    # bbk[kkk,]=as.vector(optim(rep(1,nz),Q2,lower=rep(10^-10,nz),upper=rep(10,nz),method="L-BFGS-B",gr = NULL,control=list(fnscale=-1))$par)
    bbk[kkk,]=as.vector(optim(rep(1,nz),Q2,lower=rep(1e-7,nz),
                              upper=rep(1000,nz),method="L-BFGS-B",
                              gr = NULL,control=list(fnscale=-1))$par)
    #convergence criterion  
    betat=betak[kkk,]-betak[kkk-1,]
    bt=bbk[kkk,]-bbk[kkk-1,]
    # if(sum(abs(betat))+sum(abs(bt))<10^(-3)){break}
    if(sum(abs(betat)) + sum(abs(bt)) < thresh) break
    
    #### should this be message instead?
    # print(sum(abs(betat))+sum(abs(bt)))
    
  }
  lst=list(b=bbk,beta=betak,k=kkk,a=a,g=g)
  
  # return 
  lst
}

# test
# G <- GBM_data2$G
# Y <- GBM_data2$Y
# C <- GBM_data2$C
# R2 <- GBM2_EMVS_res$R2
# Zmatrix <- Zmat_builder(R2, G)
# a0 <- 0.1
# gstr <- "scale"
# mpmath <- setup_mpmath()
# 
# NEG_em(Y, G, C, a0, gstr, Zmatrix,
#        I=10, thresh=0.001, .mpmath=mpmath)



#' Setup helper to install mpmath package
#' 
#'
#' requires reticulate
#' @export
setup_mpmath <- function(){
  # if(!reticulate::py_available()) stop("EMMultiOmics::lpcf() requires python 3.x. Please install python.")
  if(!reticulate::virtualenv_exists("emmultiomics")) reticulate::virtualenv_create("emmultiomics", packages=NULL)
  reticulate::use_virtualenv("emmultiomics")
  if(!reticulate::py_module_available("mpmath")) reticulate::virtualenv_install("emmultiomics", packages="mpmath")
  
  # return mpmath module
  reticulate::import("mpmath")
}


#' LPCF LOG Parabolic Cylinder Function
#' 
#' 
#' The function has dependency on python pacakge: mpmath
#' 
#' @param k details...
#' @param v details...
#' @param z details...
#' @param .mpmath 
#' 
#' @export
lpcf <- function(k, v, z, .mpmath){
  vz <- as.matrix(t(c(k, v, z)))
  D <- vz
  i <- D[1]
  v <- D[2:(i+1)]
  z <- D[(i+2):(2*i+1)]
  A <- rep(0, i)
  # setup mpmath for use
  # .mpmath <- .setup_mpmath()
  A <- mapply(function(x, y){
    mpmath_obj <- .mpmath$log(.mpmath$pcfd(x, y))
    as.numeric(as.character(mpmath_obj))
  }, v, z, SIMPLIFY = TRUE)
  
  # return
  A
}

#' Z-matrix builder
#' 
#' @description
#' Builds Z matrix to be used with NEG_em
#' 
#' @param R2 R-squared from EMVS result
#' @param G gene expression level
#' 
#' @examples
#' G <- GBM_data2$G
#' R2 <- GBM2_EMVS_res$R2
#' 
#' Zmat_builder(R2, G)
#' 
#' @export
Zmat_builder <- function(R2, G){
  K <- ncol(G)
  Zmatrix <- cbind(matrix(1,nrow=K,ncol=1), matrix(0,nrow=K,ncol=3))
  
  for(ww in 1:K){
    if(R2[ww] >= 0.8) Zmatrix[ww, 2] <- 1
    else{
      if(R2[ww] >= 0.2) Zmatrix[ww, 3] <- 1
      else Zmatrix[ww, 4] <- 1
    }
  }
  
  # return
  Zmatrix
}



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
#' @export
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


#' get rsq, private
.rsq <- function(x, y) summary(lm(y~x))$r.squared

#' get cindx, private
.cindx <- function(pred,actual){
  phi = 0
  phi_pred = 0
  n_input = length(actual)
  for (ci in 1:n_input){
    for (cj in 1:n_input){
      if (actual[ci]>actual[cj]){
        phi=phi+1
        if (pred[ci]>pred[cj]){
          phi_pred=phi_pred+1
        }
      }
    }
  }
  return(phi_pred/phi)
}


rm(list=ls())
# library(gridExtra)
# library(ggplot2)
# library(glmnet)
# library(MASS)
# library(cBioPortalData)
# library(dplyr)
# library(survival)
# library(coda)
# library(caret)
set.seed(2021)



### functions
# copy as local
`%dopar%` <- foreach::`%dopar%`
`%do%` <- foreach::`%do%`
`%<-%` <- future::`%<-%`

# EMVS function
# what does it return? what is the result? (read Haos paper)
# M <- what is this matrix for? data matrix? Methylation matrix
# G <- what is this G?
# grouping <- ?
# nu0? nu1? nu? lambda? a? b?
# can I be a parameter too? choose different values of I?

# clusters for parallel computing
n_cores <- parallel::detectCores() - 1
my_cluster <- parallel::makeCluster(
  n_cores,
  type="PSOCK"
)

# register them to be used with %dopar%
doParallel::registerDoParallel(cl=my_cluster)
foreach::getDoParRegistered() # TRUE
foreach::getDoParWorkers() # == n_cores

# run the function
# load rda file
load("data/GBM_data2.rda")
# load GENE_GROUP2
load("data/GENE_GROUP2.rda")


# global vars?
MECHANISTIC_MOD <- "quadratic"
CLINIC_MOD <- "em"
# we are only considering "GBM_DATA2" for now which is for quadratic mechanistic mode?
# gbm 3 and gbm 1 should be discussed later with Hao as I'm not sure what to do with them
G <- GBM_data2$G
C <- GBM_data2$C
Y <- log(GBM_data2$Y)
M <- GBM_data2$M
DELTA <- GBM_data2$Delta
grouping <- GENE_GROUP2



#### funciton
# params
nu0 <- 0.5
nu1 <- 10^3
nu <- 1
lambda <- 1
a <- 1
b <- 1


#get dimension of data
N=nrow(M)
K=ncol(G)
J=ncol(M)
#get functional cluster information
R = length(unique(grouping))
Kr = list()
for (r in unique(grouping)){
  Kr[[r]]=which(grouping==r)
}
#initialize parameters
estalpha = numeric(R)
esttau2 = numeric(R)
estOme=matrix(0,nrow=J,ncol=K)
estsig2=numeric(K)
esttheta=numeric(K)
iteration=numeric(K)

# #### export all the variables to clusters
# parallel::clusterExport(
#   cl=my_cluster,
#   varlist=c("N", "J", "K", "R", "Kr", "estalpha", "esttau2", "estOme", "estsig2", "esttheta", "iteration")
#   )
`%dopar%` <- foreach::`%dopar%`
foreach::foreach(r=1:R) %dopar%{
  # print(r)
  cat(r)
  # I=100
  I <- 10
  #initialize estimation
  ome = thresh_ome = array(0,dim=c(I,J,length(Kr[[r]])))
  sig2 = array(0,dim=c(I,length(Kr[[r]])))
  theta = array(0,dim=c(I,length(Kr[[r]])))
  alpha = array(0,I)
  tau2 = array(0,I)
  ome[1,,]=numeric(J)
  sig2[1,]=1
  theta[1,]=0.1
  alpha[1]=0.1
  tau2[1]=0.1
  #E-step
  for (i in 2:I){
    p=array(0,dim=c(J,length(Kr[[r]])))
    d=array(0,dim=c(J,length(Kr[[r]])))
    for (k in 1:length(Kr[[r]])){
      #initial value of parameters:
      for (j in 1:J){
        p1=dnorm(ome[(i-1),j,k],mean=0,sd=sqrt(sig2[i-1,k]*nu1))*theta[i-1,k]
        p0=dnorm(ome[(i-1),j,k],mean=0,sd=sqrt(sig2[i-1,k]*nu0))*(1-theta[i-1,k])
        p[j,k]=p1/(p1+p0)
        d[j,k]=(1-p[j,k])/nu0+p[j,k]/nu1
      }
      D=diag(d[,k])
      #M-step
      ome[i,,k]=coef(glmnet::glmnet(x=M,y=(G[,Kr[[r]][k]]-alpha[i-1]), intercept = F, standardize = F,
                            lambda=(sum(d[,k])/ncol(M))/N,alpha=0,penalty.factor = sqrt(d[,k])))[-1,1]
      # ome[i,,k]=(solve(D)-solve(D)%*%t(M)%*%solve(diag(nrow=N,ncol=N)+M%*%solve(D)%*%t(M))%*%M%*%solve(D))%*%t(M)%*%(G[,Kr[[r]][k]]-alpha[i-1])
      sig2[i,k]=(t(G[,Kr[[r]][k]]-M%*%ome[i,,k]-alpha[i-1])%*%(G[,Kr[[r]][k]]-M%*%ome[i,,k]-alpha[i-1])+t(sqrt(D)%*%ome[i,,k])%*%(sqrt(D)%*%ome[i,,k])+nu*lambda)/(N+J+nu+2)
      theta[i,k]=(sum(p[,k])+a-1)/(a+b+J-2)
    }
    #thresholding
    model=matrix(0,nrow=J,ncol=length(Kr[[r]]))
    for (k in 1:length(Kr[[r]])){
      postp=numeric(J)
      for (j in 1:J){
        p1=dnorm(ome[i,j,k],mean=0,sd=sqrt(sig2[i,k]*nu1))*theta[i,k]
        p0=dnorm(ome[i,j,k],mean=0,sd=sqrt(sig2[i,k]*nu0))*(1-theta[i,k])
        postp[j]=p1/(p1+p0)
        model[j,k]=ifelse(postp[j]>=0.5,1,0)
        thresh_ome[i,j,k]=ifelse(model[j,k]==1,ome[i,j,k],0)
      }
    }
    tau2[i]=(alpha[i]^2+1)/4
    alpha[i]=((tau2[i]/sig2[i,])%*%colSums(G[,Kr[[r]]]-M%*%thresh_ome[i,,]))/(1+N*sum((tau2[i]/sig2[i,])))
    # print('alpha')
    # print(alpha[i])
    cat(alpha[i])
    # print('ome:')
    # print(sum(abs(thresh_ome[i,,])))
    # print('----------------')
    if (sum(abs(ome[i,,]-ome[(i-1),,]))+sum(abs(sig2[i,]-sig2[i-1,]))+sum(abs(theta[i,]-theta[i-1,]))+abs(alpha[i]-alpha[i-1])+abs(tau2[i]-tau2[i-1])<0.0001) {
      #print('converge')
      estOme[,Kr[[r]]]=thresh_ome[i,,] 
      estsig2[Kr[[r]]]=sig2[i,] 
      esttheta[Kr[[r]]]=theta[i,] 
      esttau2[r]=tau2[i]
      estalpha[r]=alpha[i]
      iteration[Kr[[r]]]=i
      #print(paste0('alpha',estalpha[r]))
      break}
  }
  estOme[,Kr[[r]]]=thresh_ome[i,,] 
  estsig2[Kr[[r]]]=sig2[i,] 
  esttheta[Kr[[r]]]=theta[i,] 
  esttau2[r]=tau2[i]
  estalpha[r]=alpha[i]
  iteration[Kr[[r]]]=i
}
#estimate random intercept of each functional cluster
#compute R2 for each gene site
R2=numeric(K)
for (kk in 1:K){
  ind=(estOme[,kk]!=0)
  R2[kk]=ifelse(sum(ind)==0,0,(summary(lm((G[,kk]-estalpha[grouping[kk]])~M[,ind]-1))$r.squared))
}
my_result <- list(estOme=estOme,estsig2=estsig2,esttheta=esttheta,
            esttau2=esttau2,estalpha=estalpha,
            iteration=iteration,R2=R2)

# stop cluster
# parallel::stopCluster(cl=my_cluster)

# 
# EMVS = function(M,G,grouping,nu0=0.5,nu1=10^3,nu=1,lambda=1,a=1,b=1){
#   
# }




###########################################################
load("data/GBM_data2_EMVS_res.Rdata")

# What is the expected output format?
EMVS_res$estOme
EMVS_res$estsig2
EMVS_res$esttheta
EMVS_res$esttau2
EMVS_res$estalpha
EMVS_res$iteration
EMVS_res$R2

##########



EMVS_res = EMVS(M,G,grouping=GENE_GROUP2,nu0=0.5,nu1=10^3,nu=1,lambda=1,a=1,b=1)
save(EMVS_res,file=paste0('data/GBM_data2_EMVS_res.Rdata'))
load("data/GBM_data2_EMVS_res.Rdata")
N=nrow(G)
L=ncol(C)
J=ncol(M)
K=ncol(G)
if (length(args)==0){
  a0 = .1
  g0 = gstr = 1/(N^2)
}else{
  a0 = as.numeric(args[3])
  gstr = args[4]
}
mypath = paste0(getwd(),'/')
if (!file.exists(mypath)){
  dir.create(file.path(mypath))
}
# plot(EMVS_res$R2)
################################################################################
#second-stage

Zmatrix=cbind(matrix(1,nrow=K,ncol=1),matrix(0,nrow=K,ncol=3))
for(ww in 1:K){
  if(EMVS_res$R2[ww]>=0.8){Zmatrix[ww,2]=1}
  else{ if(EMVS_res$R2[ww]>=0.2){Zmatrix[ww,3]=1} else{Zmatrix[ww,4]=1}}
}

###  NEG   ###
###  NEG Group ###
# colnames(G)[abs(NEG_res$beta[NEG_res$k,])>1e-5]
# NEG_res$beta[NEG_res$k,abs(NEG_res$beta[NEG_res$k,])>1e-5]

#E-M loop
#hyperparameter

n_fold =10
set.seed(123)
folds <- createFolds(1:N, k = n_fold)
final_list = list()
print("a0")
print(a0)
print('g')
print(gstr)
for (jjj in 1:n_fold){
  test_indx = folds[[jjj]]
  train_delta = Delta[-test_indx]
  test_delta = Delta[test_indx]
  if (clinic_mod=='mcmc'){
    lst = NEG_mcmc(Y=Y[-test_indx,],G=as.matrix(G[-test_indx,]),
                   C=as.matrix(C[-test_indx,]),
                   a0=a0,gstr=gstr,Zmatrix=Zmatrix,mypath=mypath)
  }else if (clinic_mod=='em'){
    lst = NEG_em(Y=Y[-test_indx],G=as.matrix(G[-test_indx,]),
                 C=as.matrix(C[-test_indx,]),
                 a0=a0,gstr=gstr,Zmatrix=Zmatrix,mypath=mypath)
  }else if (clinic_mod=='censor'){
    lst = NEG_censor(Y=train[,1],G=as.matrix(train[,2:(K0+1)]),
                     C=as.matrix(train[,(K0+2):(K+1)]),a0=a0,gstr=gstr,
                     Zmatrix=Zmatrix,mypath=mypath,Delta=train_delta)
  }
  #final_list[[indx_g]][[indx_a]][[jjj]]=lst
  final_list[[jjj]]=lst
  save(final_list,file=paste0(mypath,"10fold_",a0,"_",gstr,"_",clinic_mod,"_censor.RData"))
}
foreach(a0=c(0.1,1,10,50)) %:% 
  foreach(gstr=c('scale')) %:%
    foreach(jjj=1:n_fold) %do%{
      test_indx = folds[[jjj]]
      train_delta = Delta[-test_indx]
      test_delta = Delta[test_indx]
      if (clinic_mod=='mcmc'){
        lst = NEG_mcmc(Y=Y[-test_indx,],G=as.matrix(G[-test_indx,]),
                       C=as.matrix(C[-test_indx,]),
                       a0=a0,gstr=gstr,Zmatrix=Zmatrix,mypath=mypath)
      }else if (clinic_mod=='em'){
        lst = NEG_em(Y=Y[-test_indx],G=as.matrix(G[-test_indx,]),
                     C=as.matrix(C[-test_indx,]),
                     a0=a0,gstr=gstr,Zmatrix=Zmatrix,mypath=mypath)
      }else if (clinic_mod=='censor'){
        lst = NEG_censor(Y=train[,1],G=as.matrix(train[,2:(K0+1)]),
                         C=as.matrix(train[,(K0+2):(K+1)]),a0=a0,gstr=gstr,
                         Zmatrix=Zmatrix,mypath=mypath,Delta=train_delta)
      }
      #final_list[[indx_g]][[indx_a]][[jjj]]=lst
      final_list[[jjj]]=lst
      save(final_list,file=paste0(mypath,"10fold_",a0,"_",gstr,"_",clinic_mod,".RData"))
    }

# #load(paste0(mypath,"10fold_",a0,"_",gstr,"_em.RData"))
# # load("~/Dropbox/harvard/emvs/code/results/mcmc_dist/10fold_0.5_1.RData")
# 
# rsq <- function(x, y) summary(lm(y~x))$r.squared
# res_table = data.frame(matrix(0,n_fold,19))
# colnames(res_table) = c('size support','r2_train','r2_test',
#                         'cindex_train','cindex_test','mse_train','mse_test',
#                         'beta_AGE_lower','beta_AGE_mean','beta_AGE_upper',
#                         'beta_PRIOR_GLIOMA_lower','beta_PRIOR_GLIOMA_mean','beta_PRIOR_GLIOMA_upper',
#                         'beta_SEX_lower','beta_SEX_mean','beta_SEX_upper',
#                         'beta_PRETREATMENT_HISTORY_lower','beta_PRETREATMENT_HISTORY_mean','beta_PRETREATMENT_HISTORY_upper')
# for (jjj in 1:n_fold){
#   test_indx = folds[[jjj]]
#   output = final_list[[jjj]]
#   output$beta = data.frame(output$beta)
#   colnames(output$beta) = colnames(G)
#   estbeta = output$beta[output$k,]
#   estbeta[(K0+1):K] =colMeans(output$beta[round(1/5*output$k):output$k,(K0+1):K])
#   pred_train = as.matrix(GenDat)[-test_indx,-1]%*%as.numeric(estbeta)
#   pred_test = as.matrix(GenDat)[test_indx,-1]%*%as.numeric(estbeta)
#   res_table[jjj,'size support']=sum(abs(estbeta)>1e-5)
#   res_table[jjj,'r2_train']=rsq(pred_train,Y[-test_indx])
#   res_table[jjj,'r2_test']=rsq(pred_test,Y[test_indx])
#   res_uotable[jjj,'cindex_train']=cindx(pred_train,Y[-test_indx])
#   res_table[jjj,'cindex_test']=cindx(pred_test,Y[test_indx])
#   res_table[jjj,'mse_train']=mean((pred_train-Y[-test_indx])^2)
#   res_table[jjj,'mse_test']=mean((pred_test-Y[test_indx])^2)
#   beta.mcmc = as.mcmc(output$beta[round(1/5*output$k):output$k,(K0+1):K])
# 
#   # find HPD interval using the CODA function
#   for (clifac in c("AGE","PRIOR_GLIOMA","SEX","PRETREATMENT_HISTORY")){
#     res_table[jjj,c(paste0('beta_',clifac,'_lower'),paste0('beta_',clifac,'_upper'))]=HPDinterval(beta.mcmc)[clifac,]
#     res_table[jjj,paste0('beta_',clifac,'_mean')]=estbeta[clifac]
#   }
# }
# #
#summarize cross-validation results
final_table = data.frame(matrix(0,nrow=8,ncol=13))
box_tables = data.frame()
colnames(final_table)[1:2] =c('a','g')
colnames(final_table)[3:13] = c('size support',
                                'r2_train',
                                'r2_test',
                                'cindex_train',
                                'cindex_test',
                                'mse_train',
                                'mse_test',
                                'AGE',
                                'PRIOR_GLIOMA',
                                'SEX',
                                'PRETREATMENT_HISTORY')

ppp=1
for (gstr in c('scale',1)){
  for (a0 in c(0.1,1,10,50)){
    fpath = paste0(mypath,"10fold_",
                   a0,"_",gstr,"_",clinic_mod,".RData")
    load(fpath)
    print(length(final_list))
    res_table = data.frame(matrix(0,n_fold,11))
    colnames(res_table) = c('size support','r2_train','r2_test',
                            'cindex_train','cindex_test','mse_train','mse_test',
                            'AGE','PRIOR_GLIOMA','SEX','PRETREATMENT_HISTORY')
    selected_biomarkers = c()
    for (jjj in 1:n_fold){
      test_indx = folds[[jjj]]
      output = final_list[[jjj]]
      output$beta = data.frame(output$beta)
      colnames(output$beta) = colnames(G)
      estbeta = output$beta[output$k,]
      selected_biomarkers = unique(c(selected_biomarkers,names(estbeta)[abs(estbeta)>1e-5]))
      pred_train = cbind(G[-test_indx,],C[-test_indx,])%*%as.numeric(estbeta)
      pred_test = cbind(G[test_indx,],C[test_indx,])%*%as.numeric(estbeta)
      res_table[jjj,'size support']=sum(abs(estbeta)>1e-5)
      res_table[jjj,'r2_train']=rsq(pred_train,Y[-test_indx])
      res_table[jjj,'r2_test']=rsq(pred_test,Y[test_indx])
      res_table[jjj,'cindex_train']=cindx(pred_train,Y[-test_indx])
      res_table[jjj,'cindex_test']=cindx(pred_test,Y[test_indx])
      res_table[jjj,'mse_train']=mean((pred_train-Y[-test_indx])^2)
      res_table[jjj,'mse_test']=mean((pred_test-Y[test_indx])^2)
      res_table[jjj,c("AGE","PRIOR_GLIOMA","SEX","PRETREATMENT_HISTORY")]=estbeta[1,1001:1004]
    }
    final_table[ppp,1]=a0
    final_table[ppp,2]=ifelse(gstr==1/(N^2),'scale',gstr)
    final_table[ppp,3:13]=round(colMeans(res_table),3)
    box_table=data.frame(a=a0,g=ifelse(gstr==1/(N^2),'scale',gstr),
                         AGE=res_table$AGE,
                         PRIOR_GLIOMA=res_table$PRIOR_GLIOMA,
                         SEX=res_table$SEX,
                         PRETREATMENT_HISTORY=res_table$PRETREATMENT_HISTORY)
    box_tables = rbind(box_tables,box_table)
    ppp=ppp+1
  }
}
box_tables$group = paste0('g=',box_tables$g,',a=',box_tables$a)
g1=ggplot(data=box_tables[1:50,],aes(x=as.factor(group),y=AGE))+
  geom_boxplot()+xlab(NULL)+ylab('beta age')+theme(axis.text.x = element_text(angle = 15))
g2=ggplot(data=box_tables[1:50,],aes(x=as.factor(group),y=PRIOR_GLIOMA))+
  geom_boxplot()+xlab(NULL)+ylab('beta prior glioma')+theme(axis.text.x = element_text(angle = 15))
g3=ggplot(data=box_tables[1:50,],aes(x=as.factor(group),y=SEX))+
  geom_boxplot()+xlab(NULL)+ylab('beta sex')+theme(axis.text.x = element_text(angle = 15))
g4=ggplot(data=box_tables[1:50,],aes(x=as.factor(group),y=PRETREATMENT_HISTORY))+
  geom_boxplot()+xlab(NULL)+ylab('beta pretreatment history')+theme(axis.text.x = element_text(angle = 15))
grid.arrange(g1,g2,g3,g4, ncol = 2)
#
final_table
paste0(final_table[1,1:7],collapse = '&')

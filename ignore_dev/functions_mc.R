###### setup
# get data
# load rda file
load("data/GBM_data2.rda")
# load GENE_GROUP2
load("data/GENE_GROUP2.rda")
# G <- GBM_data2$G
# C <- GBM_data2$C
# Y <- log(GBM_data2$Y)
# M <- GBM_data2$M
# DELTA <- GBM_data2$Delta
# grouping <- GENE_GROUP2
# EMVS ref from Hao, for ref
# load("data/GBM_data2_EMVS_res.Rdata")
# EMVS_res$R2 |> plot()

########## DEV ################## 

### set up for NEG
load("data/GBM_data2.rda")
load("data/GBM2_EMVS_res.rda")
load("data/GENE_GROUP2.rda")
# should this part be functionalized?
# name for the function? (ask Hao?)
G <- GBM_data2$G
Y <- GBM_data2$Y
C <- GBM_data2$C

N <- nrow(G)
K <- ncol(G)
L <- ncol(C)
J <- ncol(GBM_data2$M)
DELTA <- GBM_data2$Delta
R2 <- GBM2_EMVS_res$R2

# initial params
a0 <- 0.1
g0 <- gstr <- 1/N^2
# mypath can be a user argument. Tell them that destination should be a txt file
mypath <- paste0(getwd(),"/")

# construct Z matrix
Zmatrix <- cbind(matrix(1,nrow=K,ncol=1), matrix(0,nrow=K,ncol=3))
# EMVS_res R2, should be specified. it is a vector
# nrow(Zmatrix) should == length(R2) else stop
for(ww in 1:K){
  if(R2[ww] >= 0.8) Zmatrix[ww, 2] <- 1
  else{
    if(R2[ww] >= 0.2) Zmatrix[ww, 3] <- 1
    else Zmatrix[ww, 4] <- 1
  }
}

# E-M loop
# hyper params
n_fold <- 10
# set.seed(123) # let the user choose
folds <- caret::createFolds(1:N, k=n_fold)
folds
final_list <- list()
print("a0")
print(a0)
print("g")
print(g0)
print(gstr)

test_indx <- folds[[1]]
train_delta <- DELTA[-test_indx]
test_delta <- DELTA[test_indx]



#### run NEG_em
a <- a0
# outputs empty thing
write.table(" ",
            paste0(mypath,"out1",a,"_",gstr,".txt"),
            append=FALSE,row.names=FALSE,col.names=FALSE) # <- this can be converted into a table or matrix
# params
alpha <- 1
gam <- 1
s1 <- 1
s2 <- 1
.c <- 1
d <- 1
N <- nrow(GBM_data2$G)
K0 <- ncol(GBM_data2$G)
L <- ncol(GBM_data2$C)
K <- K0 + L
I <- 300
# Store the values of expectations in E-step
Elambdaj <- as.vector(rep(0,K0))
Elambdaj2 <- as.vector(rep(0,K0)) 
nz <- length(Zmatrix[1,])
# initial values
sigmak=rep(1,I)
betak=matrix(0,nrow=I,ncol=K)
bbk=matrix(1,nrow=I,ncol=nz)
v1=as.vector(rep((-(2*a+1+s1)),K0))
v2=as.vector(rep((-(2*a+1+s2)),K0))
vb=as.vector(rep(-2*(a+0.5),K0))

if(gstr=="scale"){
  g0 <- 1/N^2
}else{
  g0 <- as.numeric(gstr)
}

g <- g0

#starting iteration
# assuming iter only first
kkk <- 2
print(bbk[kkk-1,])
hf <- 1/((Zmatrix)%*%bbk[kkk-1,])
ZZ <- as.vector(abs(betak[kkk-1,1:K0])*sqrt(hf)/(sqrt(g)*sigmak[kkk-1]))

### now lpcf
# lPCFb=lpcf(K0,vb,ZZ)
vz <- as.matrix(t(c(K0, vb, ZZ))) # can just literally use this
# output
# write.table(vz,paste0(mypath,"zv",a,"_",gstr,".txt"),sep=',',append=FALSE,row.names=FALSE,col.names=FALSE) # <- this can be converted into a table or matrix
# system(paste0("python3 ", mypath, "lpcf.py", " ", mypath, " ", a," ",gstr)) <- hwat is this doing?
command_str <- paste0("python3 ", mypath, "python/lpcf.py", " ", mypath, " ", a," ",gstr)
command_str
# sys.argv[0] is the path C:/Users/mcho1/Desktop/umn/R/EMMultiOmics/python/lpcf.py
# sys.argv[1] is "mypath" in lpcf.py which is "C:/Users/mcho1/Desktop/umn/R/EMMultiOmics/" same as getwd() in R
# sys.argv[2] is a0 value, which is 0.1 (for the first iter)
# sys.argv[3] is gstr which is 4.16233090530697e-05 (for the first iter)
# D = pd.read_csv("".join(["zv",str(a0),"_",gstr,".txt"])) <- this can be replaced by vz which we already constructed earlier
D <- vz
# D |> View()
# i = int(D.columns.tolist()[0]) in lpcf.py would be in R as below
i <- D[1]
####### this we can skip maybe?
# result = []
# with open(file=zv_path, mode="r") as f:
#   for line in f:
#     result.append(list(map(float, line.split(","))))
# 
############

# v = result[0][1:(i+1)]
v <- D[2:(i+1)]
# D[2:(i+1)] |> length()
z <- D[(i+2):(2*i+1)]
# D[(i+2):(2*i+1)] |> length()
A <- rep(0, i)


############ now everything is setup for mpmath.pcfd()
# import mpmath using reticulate
# create a new environment 
# do this only not existing
if(!reticulate::virtualenv_exists("emmultiomics")) reticulate::virtualenv_create("emmultiomics")
# install mpmath if not available
# reticulate::py_module_available("mpmath")
reticulate::virtualenv_install("emmultiomics", "mpmath")
reticulate::use_virtualenv("emmultiomics")
# import
# main <- reticulate::import_main()
mpmath <- reticulate::import("mpmath")

#### this setup should be done only once.

# mpmath$pcfd(v[1], z[1]) |> as.character() |> as.numeric() |> log()
# loop for v and z
A <- mapply(function(x, y){
  mpmath_obj <- mpmath$log(mpmath$pcfd(x, y))
  as.numeric(as.character(mpmath_obj))
}, v, z, SIMPLIFY = TRUE)


A
##### so now we have the A we want, biggest challenge would be to getting this reticulate thing imported on
### the user's end
### is there a way to set this up as part of the installation for package?
### or is it possible to simply copy this generated function and save it as rda for use?
### TRY the rda method?
# foo1 <- mpmath$pcfd
# foo2 <- mpmath$log
# foo1(v[1], z[1])
# save(foo1, file="foo1.rda")
#### saving as rda doesnt work
#### it is a pointer to the python function. doesnt necessarily copy the code
#### so a few things - need to set up reticulate, the virtual env, and mpmath correctly
#### then we can use it

reticulate::virtualenv_list()
reticulate::virtualenv_remove(envname = "emmultiomics" 
                              # packages = NULL, 
                              # confirm = interactive()
                              )
#### functionalize lpcf.py in r
# set up mpmath for user
setup_mpmath <- function(){
  if(!reticulate::virtualenv_exists("emmultiomics")) reticulate::virtualenv_create("emmultiomics", packages=NULL)
  reticulate::use_virtualenv("emmultiomics")
  if(!reticulate::py_module_available("mpmath")) reticulate::virtualenv_install("emmultiomics", packages="mpmath")
  
  # return mpmath module
  reticulate::import("mpmath")
}

# test
mpmath <- setup_mpmath()
mpmath$pcfd
mpmath$log

# functionalize lpcf.py
lpcf <- function(k, v, z){
  vz <- as.matrix(t(c(k, v, z)))
  D <- vz
  i <- D[1]
  v <- D[2:(i+1)]
  z <- D[(i+2):(2*i+1)]
  A <- rep(0, i)
  # setup mpmath for use
  .mpmath <- setup_mpmath()
  A <- mapply(function(x, y){
    mpmath_obj <- .mpmath$log(.mpmath$pcfd(x, y))
    as.numeric(as.character(mpmath_obj))
  }, v, z, SIMPLIFY = TRUE)
  
  # return
  A
}

# test
lPCFb <- lpcf(K0, vb, ZZ)
lPCFb
lPCF1 <- lpcf(K0, v1, ZZ)
lPCF1
lPCF2 <- lpcf(K0, v2, ZZ)
lPCF2

# resume with NEG_em steps. ##### PICK UP FROM HERE
# set I
I <- 10

for(kkk in 2:I){
  #E-step
  print(bbk[kkk-1,])
  hf=1/((Zmatrix)%*%bbk[kkk-1,])
  ZZ=as.vector(abs(betak[kkk-1,1:K0])*sqrt(hf)/(sqrt(g)*sigmak[kkk-1]))
  lPCFb=lpcf(K0,vb,ZZ)
  lPCF1=lpcf(K0,v1,ZZ)
  lPCF2=lpcf(K0,v2,ZZ)
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
                        intercept=FALSE, standardize=F, alpha=1,
                        family="gaussian", lambda=2*sqrt(2)/sqrt(g)*sigmak[kkk-1]*(sum(Elambdaj[1:K0])/N)/(2*N),
                        penalty.factor =as.vector(Elambdaj[1:K0]))
  betak[kkk,1:K0]=coef(Bayelsso_gen)[-1]
  # Bayelsso_gen=lqa.default(x=as.matrix(G),
  #                          y=Y-as.vector(as.matrix(C)%*%betak[kkk-1,(K0+1):K]),
  #                          intercept=FALSE,
  #                          family=gaussian(),
  #                          penalty=adaptive.lasso(lambda=2*sqrt(2)/sqrt(g)*sigmak[kkk-1],
  #                                                 al.weights=as.vector(Elambdaj[1:K0])))
  # betak[kkk,1:K0]=Bayelsso_gen$coefficients
  Bayelsso_clinical = glmnet::glmnet(x=as.matrix(C),
                             y=Y-as.vector(as.matrix(G)%*%betak[kkk,1:K0]),
                             intercept=FALSE, standardize=F, alpha=0,
                             family="gaussian", lambda=1/N)
  betak[kkk,(K0+1):K]=coef(Bayelsso_clinical)[-1]
  
  plot(betak[kkk,],main=kkk)
  #update sigma (correct the typo on denominator in the paper)
  sigmak[kkk]=as.numeric(((sqrt(2)*sum(abs(Elambdaj[1:K0]/sqrt(g)*betak[kkk,1:K0]))+sqrt(2*(sum(abs(Elambdaj[1:K0]/sqrt(g)*betak[kkk,1:K0])))^2+4*((t(Y-cbind(G,C)%*%betak[kkk,])%*%(Y-cbind(G,C)%*%betak[kkk,]))^1+betak[kkk,(K0:K)]%*%betak[kkk,(K0:K)]+2*d)*(N+K+2*.c+2+L)))/(2*(N+K+2*.c+2+L))))
  #sigmak[kkk]=as.numeric(((sqrt(2)*sum(abs(Elambdaj[1:K0]/sqrt(g)*betak[kkk,1:K0]))+sqrt(2*(sum(abs(Elambdaj[1:K0]/sqrt(g)*betak[kkk,1:K0])))^2+4*((t(Y_full[train_indx]-as.matrix(train[,-1])%*%betak[kkk,])%*%(Y_full[train_indx]-as.matrix(train[,-1])%*%betak[kkk,]))^1+betak[kkk,(K0:K)]%*%betak[kkk,(K0:K)]+2*d)*(N+K+2*c+2+L)))/(2*(N+K+2*c+2+L))))
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
  if(sum(abs(betat))+sum(abs(bt))<10^(-3)){break}
  print(sum(abs(betat))+sum(abs(bt)))
}


## lpcf.py



# lst <- NEG_em(Y=Y[-test_indx], G=as.matrix(G[-test_indx, ]),
#               C=as.matrix(C[-test_indx]),
#               a0=a0, gstr=gstr, Zmatrix=Zmatrix, mypath=mypath)
# final_list <- lst

# 
# for(jjj in 1:n_fold){
#   test_indx <- folds[[jjj]]
#   train_delta <- DELTA[-test_indx]
#   test_delta <- DELTA[test_indx]
#   # assumme clinic mod = "em" only
#   lst <- NEG_em(Y=Y[-test_indx], G=as.matrix(G[-test_indx, ]),
#                 C=as.matrix(C[-test_indx]),
#                 a0=a0, gstr=gstr, Zmatrix=Zmatrix, mypath=mypath)
#   final_list <- lst
#   save(final_list,file=paste0(mypath,"10fold_",a0,"_",gstr,"_",clinic_mod,"_censor.RData"))
# }

##### 
second_stage <- function(){
  NULL
}


# second stage 

#### NEG args
# assume default values are from when there is no commandargs
####### funcitonalize probably below?
a0 <- 0.1
N <- nrow(G)
K <- ncol(G)
L <- ncol(C)
J <- ncol(M)

g0 <- 1/(N^2)
gstr <- 1/(N^2) # used to name saved output files. do we need to keep?

Zmatrix <- cbind(matrix(1, nrow=K, ncol=1), matrix(0, nrow=K, ncol=3))
for(ww in 1:K){
  if(EMVS_res$R2[ww] >= 0.2) Zmatrix[ww, 2] <- 1
  else{
    if(EMVS_res$R2[ww]>=0.2) Zmatrix[ww,3] <- 1 
    else Zmatrix[ww,4] <- 1
  }
}
n_fold <- 10
caret::createFolds(1:N, k=n_fold)
final_list <- list()

for(jjj in 1:n_fold){
  test_indx <- folds[[jjj]]
  train_delta <- DELTA[-test_indx]
  test_delta <- DELTA[test_indx]
  lst = NEG_em(Y=Y[-test_indx],G=as.matrix(G[-test_indx,]),
               C=as.matrix(C[-test_indx,]),
               a0=a0,gstr=gstr,Zmatrix=Zmatrix,mypath=mypath)
  final_list[[jjj]]=lst
  
}
####### funcitonalize probably above






##################################################### DEV for second stage ##########################
# check if the user has set up mpmath? or NEG_em is going to simply call this object by default?
.mpmath <- setup_mpmath()


### set up for second stage modeling
load("data/GBM_data2.rda")
load("data/GBM2_EMVS_res.rda")
# load("data/GENE_GROUP2.rda")



# construct Z matrix
# K is from K=ncol(G), G has to be parameterized then?, user param
K <- ncol(GBM_data2$G)
# R2 is from EMVS result, user parameter
R2 <- GBM2_EMVS_res$R2

Zmatrix <- cbind(matrix(1,nrow=K,ncol=1), matrix(0,nrow=K,ncol=3))

for(ww in 1:K){
  if(R2[ww] >= 0.8) Zmatrix[ww, 2] <- 1
  else{
    if(R2[ww] >= 0.2) Zmatrix[ww, 3] <- 1
    else Zmatrix[ww, 4] <- 1
  }
}

# E-M loop step
n_fold <- 10 # user parameter
set.seed(123) # user parameter
# N is from nrow(G)
N <- nrow(GBM_data2$G) # user param should be G
folds <- caret::createFolds(1:N, k=n_fold)
final_list <- list() # to store result
Delta <- GBM_data2$Delta # user parameter

# at the first loop Hao only used a0 = 0.1 but later he used c(0.1, 1, 10, 50) so this second stage function
# should be loopable with diff params of a0 and gstr
a0 <- 0.1 # user param
gstr <- 1/N^2 # user param <- from the txt files Hao created gstr values he used were c(1, 4.162331e-05, "scale") which is the same as c(1, 1/N^2, 1/(N^2)) so scale i just the same as 1/N^2
# but not sure where to find these from the code he wrote main.R
Y <- GBM_data2$Y # user param
G <- GBM_data2$G # user param
C <- GBM_data2$C # user param

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
                I=10,
                thresh=0.001,
                .mpmath=.mpmath)
  # save result
  final_list[[jjj]] <- lst
}



##### Built 2nd stage helepr, now need to create the result data frames and graphs
# set up
G <- GBM_data2$G
Y <- GBM_data2$Y
C <- GBM_data2$C
Delta <- GBM_data2$Delta
R2 <- GBM2_EMVS_res$R2
a0 <- 0.1
gstr <- "scale"
n_fold <- 10
random_seed <- 123

# first part
result1 <- second_stage_helper(Y, G, C, Delta, R2, a0, gstr, n_fold, random_seed,
                    I=10, thresh=0.001)

# multiple a0
a0 <- c(0.1, 1, 10, 50)
result2 <- lapply(a0, function(a) second_stage_helper(Y, G, C, Delta, R2, a0=a, gstr, n_fold, random_seed, 
                                           I=10, thresh=0.001)) |> 
  setNames(a0)



# save results temporarily and come back for dev
# save(result1, file="ignore_dev/2nd_stage_result1.rda")
# save(result2, file="ignore_dev/2nd_stage_result2.rda")
##### PICK UP FROM HERE - develop summary functions (tables and plot functions?)
# load results
load("ignore_dev/2nd_stage_result1.rda")
load("ignore_dev/2nd_stage_result2.rda")

result1 # result from 2nd stage helper
result2 # result from 2nd stage helper with loop over a0

### Hao's code
# initialize result table
cols <- c("a", "g", "size_support",
          "r2_train", "r2_test", "cindex_train", "cindex_test",
          "mse_train", "mse_test", "AGE", "PRIOR_GLIOMA",
          "SEX", "PRETREATMENT_HISTORY") # does this need to be parameterized? not sure what hes doing here
final_table <- data.frame(matrix(0, nrow=8, ncol=length(cols))) # where is the nrow and ncol from?
colnames(final_table) <- cols
final_table
box_tables <- data.frame()
box_tables


# colnames(final_table)[1:2] <- c('a','g')
final_table
# colnames(final_table)[3:13] = c('size support',
#                                 'r2_train',
#                                 'r2_test',
#                                 'cindex_train',
#                                 'cindex_test',
#                                 'mse_train',
#                                 'mse_test',
#                                 'AGE',
#                                 'PRIOR_GLIOMA',
#                                 'SEX',
#                                 'PRETREATMENT_HISTORY')

### read in 10 fold results for comparison
# load("ignore_dev/hao_outputs/10fold_0.1_1_em.RData")
# load("ignore_dev/hao_outputs/10fold_0.1_4.16233090530697e-05_em_censor.RData")
# load("ignore_dev/hao_outputs/10fold_0.1_4.16233090530697e-05_em.RData")
# print(length(final_list))

n_fold <- 10 # borrow from 2nd stage helper or make this user input?
# or potentially combine the cross validation and summary fx together
###### just use length(final_list)
# plotting can be on its own?
res_cols <- c("size_support", "r2_train", "r2_test", "cindex_train", "cindex_test",
              "mse_train", "mse_test", 
              # covariates that should be specific to the dataset? - ask Hao
              "AGE", "PRIOR_GLIOMA", "SEX", "PRETREATMENT_HISTORY")
res_table <- data.frame(matrix(0, n_fold, length(res_cols))) |> 
  setNames(res_cols)
res_table




# ### EMVS parallel backend?
# EMVS_par <- function(M,G,grouping,nu0=0.5,nu1=10^3, nu=1, lambda=1, a=1, b=1, I=100,
#                  THRESH=0.0001){
#   #get dimension of data
#   N=nrow(M)
#   K=ncol(G)
#   J=ncol(M)
#   #get functional cluster information
#   R = length(unique(grouping))
#   Kr = list()
#   for (r in unique(grouping)){
#     Kr[[r]]=which(grouping==r)
#   }
#   #initialize parameters
#   estalpha = numeric(R)
#   esttau2 = numeric(R)
#   estOme=matrix(0,nrow=J,ncol=K)
#   estsig2=numeric(K)
#   esttheta=numeric(K)
#   iteration=numeric(K)
#   
#   # parallel processing
#   n_cores <- parallel::detectCores() - 1
#   cl <- parallel::makeCluster(n_cores)
#   # `%dopar%` <- foreach::`%dopar%`
#   parallel::clusterExport(cl, envir=environment(), varlist=c(
#     "M", "G", "grouping",
#     "N", "J", "K", "R", "Kr",
#     "estalpha", "esttau2", "estOme", "estsig2", "esttheta", "iteration",
#     "nu0", "nu1", "nu", "lambda", "a", "b", "I", "THRESH"
#   ))
#   
#   result <- parallel::parLapply(cl, 1:R, function(r){
#     # track progress
#     # pb <- txtProgressBar(min=0, max=R, initial=0, style=3)
#     ome = thresh_ome = array(0,dim=c(I,J,length(Kr[[r]])))
#     sig2 = array(0,dim=c(I,length(Kr[[r]])))
#     theta = array(0,dim=c(I,length(Kr[[r]])))
#     alpha = array(0,I)
#     tau2 = array(0,I)
#     ome[1,,]=numeric(J)
#     sig2[1,]=1
#     theta[1,]=0.1
#     alpha[1]=0.1
#     tau2[1]=0.1
#     #E-step
#     for (i in 2:I){
#       p=array(0,dim=c(J,length(Kr[[r]])))
#       d=array(0,dim=c(J,length(Kr[[r]])))
#       for (k in 1:length(Kr[[r]])){
#         #initial value of parameters:
#         for (j in 1:J){
#           p1=dnorm(ome[(i-1),j,k],mean=0,sd=sqrt(sig2[i-1,k]*nu1))*theta[i-1,k]
#           p0=dnorm(ome[(i-1),j,k],mean=0,sd=sqrt(sig2[i-1,k]*nu0))*(1-theta[i-1,k])
#           p[j,k]=p1/(p1+p0)
#           d[j,k]=(1-p[j,k])/nu0+p[j,k]/nu1
#         }
#         D=diag(d[,k])
#         #M-step
#         ome[i,,k]=coef(glmnet::glmnet(x=M,y=(G[,Kr[[r]][k]]-alpha[i-1]), intercept = F, standardize = F,
#                                       lambda=(sum(d[,k])/ncol(M))/N,alpha=0,penalty.factor = sqrt(d[,k])))[-1,1]
#         # ome[i,,k]=(solve(D)-solve(D)%*%t(M)%*%solve(diag(nrow=N,ncol=N)+M%*%solve(D)%*%t(M))%*%M%*%solve(D))%*%t(M)%*%(G[,Kr[[r]][k]]-alpha[i-1])
#         sig2[i,k]=(t(G[,Kr[[r]][k]]-M%*%ome[i,,k]-alpha[i-1])%*%(G[,Kr[[r]][k]]-M%*%ome[i,,k]-alpha[i-1])+t(sqrt(D)%*%ome[i,,k])%*%(sqrt(D)%*%ome[i,,k])+nu*lambda)/(N+J+nu+2)
#         theta[i,k]=(sum(p[,k])+a-1)/(a+b+J-2)
#       }
#       #thresholding
#       model=matrix(0,nrow=J,ncol=length(Kr[[r]]))
#       for (k in 1:length(Kr[[r]])){
#         postp=numeric(J)
#         for (j in 1:J){
#           p1=dnorm(ome[i,j,k],mean=0,sd=sqrt(sig2[i,k]*nu1))*theta[i,k]
#           p0=dnorm(ome[i,j,k],mean=0,sd=sqrt(sig2[i,k]*nu0))*(1-theta[i,k])
#           postp[j]=p1/(p1+p0)
#           model[j,k]=ifelse(postp[j]>=0.5,1,0)
#           thresh_ome[i,j,k]=ifelse(model[j,k]==1,ome[i,j,k],0)
#         }
#       }
#       tau2[i]=(alpha[i]^2+1)/4
#       alpha[i]=((tau2[i]/sig2[i,])%*%colSums(G[,Kr[[r]]]-M%*%thresh_ome[i,,]))/(1+N*sum((tau2[i]/sig2[i,])))
#       # print('alpha')
#       # print(alpha[i])
#       # print('ome:')
#       # print(sum(abs(thresh_ome[i,,])))
#       # print('----------------')
#       if (sum(abs(ome[i,,]-ome[(i-1),,]))+sum(abs(sig2[i,]-sig2[i-1,]))+sum(abs(theta[i,]-theta[i-1,]))+abs(alpha[i]-alpha[i-1])+abs(tau2[i]-tau2[i-1])<THRESH) { # can this threshold value be chosen by user?
#         #print('converge')
#         estOme[,Kr[[r]]]=thresh_ome[i,,] 
#         estsig2[Kr[[r]]]=sig2[i,] 
#         esttheta[Kr[[r]]]=theta[i,] 
#         esttau2[r]=tau2[i]
#         estalpha[r]=alpha[i]
#         iteration[Kr[[r]]]=i
#         #print(paste0('alpha',estalpha[r]))
#         break}
#     }
#     estOme[,Kr[[r]]]=thresh_ome[i,,] 
#     estsig2[Kr[[r]]]=sig2[i,] 
#     esttheta[Kr[[r]]]=theta[i,] 
#     esttau2[r]=tau2[i]
#     estalpha[r]=alpha[i]
#     iteration[Kr[[r]]]=i
#     
#     # progress tracker
#     # setTxtProgressBar(pb, r)
#     
#     return(list(estOme=estOme,estsig2=estsig2,esttheta=esttheta,
#                 esttau2=esttau2,estalpha=estalpha,
#                 iteration=iteration))
#   })
#   
#   # stop cluster
#   parallel::stopCluster(cl)
#   # gcinfo(verbose=FALSE)
#   # gc(verbose=FALSE)
#   
#   
#   #estimate random intercept of each functional cluster
#   #compute R2 for each gene site
#   R2=numeric(K)
#   for (kk in 1:K){
#     ind=(estOme[,kk]!=0)
#     R2[kk]=ifelse(sum(ind)==0,0,(summary(lm((G[,kk]-estalpha[grouping[kk]])~M[,ind]-1))$r.squared))
#   }
#   
#   result$R2 <- R2
#   # return
#   result
# }
# 
# 
# EMVS_par(M, G, grouping, I=2, THRESH=0.1)
# 
# 
# 
# 
# 



########## summary and plotting
load("ignore_dev/2nd_stage_result1.rda")
load("ignore_dev/2nd_stage_result2.rda")

final_table <- data.frame(matrix(0, nrow=8, ncol=13)) # 8 and 13 coming from?
final_table
box_tables <- data.frame()
box_tables







############################################# start off fresh here ######
devtools::document()
devtools::build()
library(EMMultiOmics)

# run cv helper
G <- GBM_data2$G
Y <- GBM_data2$Y
C <- GBM_data2$C
Delta <- GBM_data2$Delta
R2 <- GBM2_EMVS_res$R2
a0 <- 0.1
gstr <- "scale"
n_fold <- 10
random_seed <- 123

NEG_list1 <- cv_helper(Y, G, C, Delta, R2, a0, gstr, n_fold, random_seed,
                       I=10, thresh=0.001)
# multiple a0
a0 <- c(0.1, 1, 10, 50)
NEG_list2 <- lapply(a0, function(a) cv_helper(Y, G, C, Delta, R2, a0=a, gstr, n_fold, random_seed, 
                                 I=10, thresh=0.001)) |> 
  setNames(a0)

#### now go back to developing cv helper?
# make folds
####
n_fold = 10
set.seed(123)
N <- nrow(G)
folds <- caret::createFolds(1:N, k = n_fold)
####


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
for(gstr in c("scale", 1)){
  for (a0 in c(0.1, 1, 10, 50)){
    res_cols <- c("size_support", "r2_train", "r2_test", "cindex_train", "cindex_test",
                  "mse_train", "mse_test", 
                  # covariates that should be specific to the dataset? - ask Hao
                  "AGE", "PRIOR_GLIOMA", "SEX", "PRETREATMENT_HISTORY")
    res_table <- data.frame(matrix(0, n_fold, length(res_cols))) |> 
      setNames(res_cols)
    
    selected_biomarkers <- c()
    
    for(jjj in 1:n_fold){
      test_indx = folds[[jjj]]
      output = NEG_list1[[jjj]]
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
    final_table[ppp,1]=a0
    final_table[ppp,2]=ifelse(gstr==1/(N^2),'scale',gstr)
    final_table[ppp,3:13]=round(colMeans(res_table),3) # 3:13 from where?
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


final_table



###### TRY
# run cv helper
G <- GBM_data2$G
Y <- GBM_data2$Y
C <- GBM_data2$C
Delta <- GBM_data2$Delta
R2 <- GBM2_EMVS_res$R2
a0 <- 0.1
gstr <- "scale"
n_fold <- 10
random_seed <- 123

final_list1 <- cv_helper(Y, G, C, Delta, R2, a0, gstr, n_fold, random_seed,
                       I=10, thresh=0.001)

# multiple a0 and gstr
a0 <- c(0.1, 1, 10, 50)
gstr <- c("scale", 1)
result2 <- lapply(
  gstr, 
  function(g) lapply(
    a0, 
    function(a) cv_helper(Y, G, C, Delta, R2, a0=a, gstr=g, n_fold, random_seed, I=10, thresh=0.001)
  ))

# summary plot
summary_plot(box_tables=final_list2$box_table)














###### 1. MultiOmics
# start with EMVS
# set up fx params
devtools::load_all()
# params for EMVS
M <- GBM_data2$M
G <- GBM_data2$G
# gene_grouping <- GENE_GROUP2
fp <- system.file("eg_fx_classification.txt", package="EMMultiOmics")
gene_grouping <- get_grouping(eg_gene_symbols, fp)
# gene_grouping
### include ... to specify
# nu0, nu1, nu, lambda, a, b
# I and thresh, make separate input for EMVS? like EMVS.I, EMVS.thresh

# params for second stage
G



# within MultiOmics, run EMVS
EMVS_result <- EMVS(M,G, grouping=gene_grouping, I=3, thresh=0.001)

### 2nd stage
# params
Y <- GBM_data2$Y
C <- GBM_data2$C
R2 <- EMVS_result$R2
Delta <- GBM_data2$Delta
.mpmath=setup_mpmath()

a0 <- 0.1
gstr <- 1
n_fold <- 10
random_seed <- 123



# run 2nd stage
N <- nrow(G)
K <- ncol(G)
J <- ncol(M)
L <- ncol(C)

# build Z matrix
Zmatrix <- Zmat_builder(R2, G)

# E-M loop step
# run cross validation
set.seed(random_seed)
folds <- caret::createFolds(1:N, k=n_fold)
# NEG result list
# NEG_list <- list() 

# initialize result objects
selected_biomarkers <- c()
res_cols <- c("n_selected_vars", "r2_train", "r2_test",
              "cindex_train", "cindex_test", "mse_train", "mse_test",
              colnames(C))

res_table <- data.frame(matrix(0, n_fold, length(res_cols))) |> 
  setNames(res_cols)

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
                       I=3, # make it part of ...? # or separate parameter, NEG.I
                       thresh=0.001, # make it part of ...? # or separate parameter, NEG.thresh
                       # ...,
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
  
}

# take the average
performance_df <- data.frame(a=a0, #g=ifelse(gstr==1/(N^2),'scale',gstr),
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

# calculate SE?
# reference https://stats.stackexchange.com/questions/44838/how-are-the-standard-errors-of-coefficients-calculated-in-a-regression/44841#44841
c_beta <- estbeta[colnames(C)] |> as.matrix() |> t() |> `colnames<-`("beta")
c_beta

sigma_sq <- sum((Y - C %*% c_beta)^2) / (nrow(C) - ncol(C))
# vcov_mat <- sigma_sq * chol2inv(chol(t(C) %*% C)) 
vcov_mat <- sigma_sq * solve(t(C) %*% C)
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

print(final_result)

final_result

# use the quantiles of beta instead?

# try whatever sounak suggested
CCt <- C %*% t(C)
# (CCt + Iq) inverse
mat1 <- solve(CCt + diag(x=1, nrow=dim(CCt)[1], ncol=dim(CCt)[2])) |>
  `colnames<-`(NULL) |>
  `rownames<-`(NULL)

mat1
sqrt(mat1 |> diag())
# mat2 <- MASS::ginv(CCt + diag(x=1, nrow=dim(CCt)[1], ncol=dim(CCt)[2]))
# mat2
mat1 %*% t(C)


# plot
# where do we get colnames C from?
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

# how final_result# how about full design matrix?
# full_beta <- estbeta[, estbeta != 0] |> as.matrix() |> t()
# nonzero_estbeta <- estbeta[, estbeta != 0]
# full_dx_mat <- cbind(G[, colnames(G) %in% names(nonzero_estbeta)], C[, colnames(C) %in% names(nonzero_estbeta)])
# full_sigma_sq <- sum((Y - full_dx_mat %*% full_beta)^2) / (nrow(full_dx_mat) - ncol(full_dx_mat))
# full_vcov_mat <- full_sigma_sq * chol2inv(chol(t(full_dx_mat) %*% full_dx_mat))
# full_std_err <- sqrt(diag(full_vcov_mat))


#### multiOmics_sensitivity
# loop over g and a values
a_values <- c()
g_values <- c()











# 
# 
# 
# ### not sure if I'm estimating this correctly
# 
# 
# # take the last output and develop
# # this is happening inside the loop
# # summarize NEG result
# output <- NEG_list[[10]]
# output$beta <- data.frame(output$beta) # nrow of this dataframe = I in NEG_em
# colnames(output$beta) <- c(colnames(G), colnames(C))
# # output$beta is all of the beta estimates
# estbeta <- output$beta[output$k, ] # take the last iterations row as the estimated beta
# 
# # selected biomarkers (or the genes)
# names(estbeta)[abs(estbeta) > 1e-5] # this has to be repeated for each a and gstr values
# selected_biomarkers <- unique(c(selected_biomarkers, names(estbeta)[abs(estbeta) > 1e-5]))
# 
# # ask this, should we drop the clinical variables? consider only the number of genes in the selected biomarkers?
# selected_biomarkers <- unique(c(selected_biomarkers, 
#                                 names(estbeta)[abs(estbeta) > 1e-5 & !names(estbeta) %in% colnames(C)]))
# 
# # evaluate cross validation result
# pred_train <- cbind(G[-test_indx,], C[-test_indx,]) %*% as.numeric(estbeta)
# pred_test <- cbind(G[test_indx,], C[test_indx,]) %*% as.numeric(estbeta)
# # number of selected biomarkers is "size support"
# size_support <- sum(abs(estbeta) > 1e-5 & !names(estbeta) %in% colnames(C)) # originally it was just sum(abs(estbeta > 1e-5))
# r2_train <- .rsq(pred_train, Y[-test_indx])
# r2_test <- .rsq(pred_test, Y[test_indx])
# cindex_train <- .cindx(pred_train, Y[-test_indx])
# cindex_test <- .cindx(pred_test, Y[test_indx])
# mse_train <- mean((pred_train - Y[-test_indx])^2)
# mse_test <- mean((pred_test - Y[test_indx])^2)
# clinical_beta <- estbeta[colnames(C)]
# 
# 
# # after the loop, remove NA
# # selected_biomarkers <- selected_biomarkers[!is.na(selected_biomarkers)]
# .rsq(pred_train, Y[-test_indx])
# lm(Y[-test_indx] ~ pred_train) |> summary()
# lm(Y[test_indx] ~ pred_test) |> summary()
# 
# # squared error
# (pred_train - Y[-test_indx])^2
# (pred_test - Y[test_indx])^2
# 
# 



#### test
# params
M <- GBM_data2$M
G <- GBM_data2$G
grouping <- GENE_GROUP2
Y <- GBM_data2$Y
C <- GBM_data2$C
a0 <- 0.1
gstr <- "scale"
Delta <- GBM_data2$Delta
n_fold <- 10
random_seed <- 123
EMVS_I <- 3
NEG_I <- 3
EMVS_thresh <- 0.0001
NEG_thresh <- 0.0001

# run
multiOmics(
  M=M,
  G=G,
  grouping=grouping,
  Y=Y,
  C=C,
  a0=a0,
  gstr=gstr,
  Delta=Delta,
  n_fold=n_fold,
  random_seed = random_seed,
  EMVS_I = EMVS_I,
  NEG_I=NEG_I,
  EMVS_thresh = EMVS_thresh,
  NEG_thresh = NEG_thresh
  )








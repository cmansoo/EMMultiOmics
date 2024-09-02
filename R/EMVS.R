#' Expectation Maximization Variable Selection Function
#' 
#' Performs EMVS algorithm for genomics data.
#' The default values for `nu0`, `nu1`, `lambda`, `a`, `b` were chosen per publication: *Veronika Ročková & Edward I. George (2013): EMVS: The EM Approach to Bayesian Variable Selection, Journal of the American Statistical Association*
#' 
#' @param M DNA Methylation matrix
#' @param G Gene expression level
#' @param grouping Gene grouping
#' @param nu0 parameter 0 for spike-and-slab Gaussian mixture prior on \eqn{\beta}
#' @param nu1 parameter 1 for spike-and-slab Gaussian mixture prior on \eqn{\beta}
#' @param lambda For the prior on \eqn{\sigma2}, an inverse gamma prior \eqn{\pi(\sigma2 | \gamma) = IG(\nu/2, \nu\lambda/2)}, <- what is lambda here?
#' @param a Parameter \eqn{a} for the beta prior \eqn{\pi(\theta) \propto \theta a-1(1-\theta)b-1}, \eqn{a, b > 0}
#' @param b Parameter \eqn{b} for the beta prior \eqn{\pi(\theta) \propto \theta a-1(1-\theta)b-1}, \eqn{a, b > 0}
#' @param I Maximum number of iterations of EMVS
#' @param thresh Threshold for convergence criterion.
#' 
#' @examples
#' 
#' # example data
#' M <- GBM_data2$M
#' G <- GBM_data2$G
#' grouping <- GENE_GROUP2
#' 
#' EMVS(M, G, grouping, I=10, thresh=0.001)
#' # EMVS_par(M, G, grouping, I=10, thresh=0.001)
#' 
#' @export
EMVS <- function(M, G, grouping, nu0=0.5, nu1=10^3, nu=1, lambda=1, a=1, b=1, I=100,
                 thresh=0.0001){
  # get dimension of data
  N=nrow(M)
  K=ncol(G)
  J=ncol(M)
  # get functional cluster information
  R = length(unique(grouping))
  Kr = list()
  for (r in unique(grouping)){
    Kr[[r]]=which(grouping==r)
  }
  # initialize parameters
  estalpha = numeric(R)
  esttau2 = numeric(R)
  estOme=matrix(0,nrow=J,ncol=K)
  estsig2=numeric(K)
  esttheta=numeric(K)
  iteration=numeric(K)
  
  message("Tracking progress...")
  pb <- txtProgressBar(min=0, max=R, initial=0, style=3)
  for(r in 1:R){
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
    # E-step
    for (i in 2:I){
      p=array(0,dim=c(J,length(Kr[[r]])))
      d=array(0,dim=c(J,length(Kr[[r]])))
      for (k in 1:length(Kr[[r]])){
        # initial value of parameters:
        for (j in 1:J){
          p1=dnorm(ome[(i-1),j,k],mean=0,sd=sqrt(sig2[i-1,k]*nu1))*theta[i-1,k]
          p0=dnorm(ome[(i-1),j,k],mean=0,sd=sqrt(sig2[i-1,k]*nu0))*(1-theta[i-1,k])
          p[j,k]=p1/(p1+p0)
          d[j,k]=(1-p[j,k])/nu0+p[j,k]/nu1
        }
        D=diag(d[,k])
        # M-step
        ome[i,,k]=coef(glmnet::glmnet(x=M,y=(G[,Kr[[r]][k]]-alpha[i-1]), intercept = F, standardize = F,
                                      lambda=(sum(d[,k])/ncol(M))/N,alpha=0,penalty.factor = sqrt(d[,k])))[-1,1]
        sig2[i,k]=(t(G[,Kr[[r]][k]]-M%*%ome[i,,k]-alpha[i-1])%*%(G[,Kr[[r]][k]]-M%*%ome[i,,k]-alpha[i-1])+t(sqrt(D)%*%ome[i,,k])%*%(sqrt(D)%*%ome[i,,k])+nu*lambda)/(N+J+nu+2)
        theta[i,k]=(sum(p[,k])+a-1)/(a+b+J-2)
      }
      # thresholding
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
      summed_q <- sum(abs(ome[i,,]-ome[(i-1),,]))+sum(abs(sig2[i,]-sig2[i-1,]))+sum(abs(theta[i,]-theta[i-1,]))+abs(alpha[i]-alpha[i-1])+abs(tau2[i]-tau2[i-1])
      
      if (summed_q < thresh) {
        # print('converge')
        estOme[,Kr[[r]]]=thresh_ome[i,,] 
        estsig2[Kr[[r]]]=sig2[i,] 
        esttheta[Kr[[r]]]=theta[i,] 
        esttau2[r]=tau2[i]
        estalpha[r]=alpha[i]
        iteration[Kr[[r]]]=i
        # print(paste0('alpha',estalpha[r]))
        break}
    }
    estOme[,Kr[[r]]]=thresh_ome[i,,] 
    estsig2[Kr[[r]]]=sig2[i,] 
    esttheta[Kr[[r]]]=theta[i,] 
    esttau2[r]=tau2[i]
    estalpha[r]=alpha[i]
    iteration[Kr[[r]]]=i
    
    # track progress
    setTxtProgressBar(pb, r)
  }
  
  # estimate random intercept of each functional cluster
  # compute R2 for each gene site
  R2=numeric(K)
  for (kk in 1:K){
    ind=(estOme[,kk]!=0)
    R2[kk]=ifelse(sum(ind)==0,0,(summary(lm((G[,kk]-estalpha[grouping[kk]])~M[,ind]-1))$r.squared))
  }
  
  message("\nComplete!")
  
  # return
  list(estOme=estOme,
       estsig2=estsig2,
       esttheta=esttheta,
       esttau2=esttau2,
       estalpha=estalpha,
       iteration=iteration,
       R2=R2)
  
}














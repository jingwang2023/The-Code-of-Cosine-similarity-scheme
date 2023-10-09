library(lsa)
library(MASS)
library(depth)
library(mvtnorm)
library(lcmix)
library(corpcor)

rm(list=ls())

# m: Size of the Reference Sample
# n: Size of the Test sample
# p: Number of the variables in a sample
# nu: Partition size of reference sample
# H: Trial value of the Upper Control Limit

###### Computation of IC-MRL for IP Design for a given (m,n,p,H)  #######


######## m: reference sample size ########
######## n: test sample size      ########
######## p: dimension             ########

#######################################################################################


####################################################################

#ref0: subsample of reference sample
#ref1: Reference sample excluding subsample for benchmarking
#test: test sample

WRS_M = function(ref0,ref1,test){
  m0=nrow(ref0) #choosing benchmark
  m1=nrow(ref1)  #remain reference
  n=nrow(test)
  
  if (m0==1){
    ref0=as.vector(ref0)
    X_b_med=ref0
  }
  else{
    X_b_med = med(ref0, method = "Spatial")$median
  }
  
  mu0=colMeans(ref0)
  v=cov.shrink(ref0)
  
  
  U_ref=sapply(1:m1,function(i){mahalanobis(ref1[i,],mu0,v)})
  V_test=sapply(1:n,function(i){mahalanobis(test[i,],mu0,v)})
  
  Z = ((sum(rank(c(U_ref,V_test))[(m1+1):(m1+n)]) - n*(m1+n+1)/2 )^2 )/ (m1*n*(m1+n+1)/12)
  
  return(Z)
  
}

WRS_A = function(ref0,ref1,test){
  m0=nrow(ref0) #choosing benchmark
  m1=nrow(ref1)  #remain reference
  n=nrow(test)
  
  if (m0==1){
    ref0=as.vector(ref0)
    X_b_med=ref0
  }
  else{
    X_b_med = med(ref0, method = "Spatial")$median
  }
  
  
  U_ref=sapply(1:m1,function(i){acos(cosine(X_b_med, ref1[i,]))})
  V_test=sapply(1:n,function(i){acos(cosine(X_b_med, test[i,]))})
  
  
  Z = ((sum(rank(c(U_ref,V_test))[(m1+1):(m1+n)]) - n*(m1+n+1)/2 )^2 )/ (m1*n*(m1+n+1)/12)
  
  return(Z)
  
}

#####################################################################################

IC.MRL2=function(m,n,p,H,del1,del2,del3,rho,lim,sim.num){
  set.seed(060819)
  cov.mat=matrix(rep(rho,p*p),p,p)+diag(rep((1-rho),p))
  RL_Array=rep(0,sim.num)
  
  for(i in 1:sim.num){
    
    r=matrix(0:0,nrow=m,ncol=6)
    r[,1]=1
    ref=rmvt(m,sigma=cov.mat,df=3)+r
    
    nu=round(m/5)
    index = c(sample(1:m,nu))
    ref0=ref[index,]
    if (nu==1) ref0=t(as.matrix(ref0))
    
    ref1=ref[-index,]
    
    plot_stat=k=0
    while(plot_stat<H && k<=lim){
      k=k+1
      r1=matrix(0:0,nrow=n,ncol=6)
      r1[,1]=del1
      r1[,2]=del2
      r1[,3]=r1[,4]=r1[,5]=r1[,6]=del3
      test=rmvt(n,sigma=cov.mat,df=3)+r1
      
      Z_M = WRS_M(ref0,ref1,test)
      Z_A = WRS_A(ref0,ref1,test)
      plot_stat=max(Z_M,Z_A)
    }
    if((i%%3000)==0) print(i)
    RL_Array[i]=k
  }
  obsARL<-sum(RL_Array)/sim.num # observed average run length
  sdN=sd(RL_Array) # observed sd of run length
  quantsN<-quantile(RL_Array,probs=c(0.05,0.25,0.5,0.75,0.95)) # Percentiles of run length distribution
  print(c(m,n,nu,rho,p,H,del1,del2,del3,obsARL,sdN,as.numeric(quantsN)))
}
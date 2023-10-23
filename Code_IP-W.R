library(lsa)
library(depth)
library(mvtnorm)


rm(list=ls())

# m: Size of the Reference Sample
# n: Size of the Test sample
# p: Number of the variables in a sample
# nu: Partition size of reference sample, mu=max
# H: Trial value of the Upper Control Limit
#ref0: subsample of reference sample
#ref1: Reference sample excluding subsample for benchmarking
#test: test sample

IPLep2 = function(ref0,ref1,test){
  m0=nrow(ref0) 
  m1=nrow(ref1)  
  n=nrow(test)
  
  U_ref=sapply(1:m1,function(i){
    c1=rbind(ref0,ref1[i,])
    dist(c1,method = "euclidean")[1]})
  
  V_test=sapply(1:n,function(i){
    c1=rbind(ref0,test[i,])
    dist(c1,method = "euclidean")[1]})
  
  
  Z = ((sum(rank(c(U_ref,V_test))[(m1+1):(m1+n)]) - n*(m1+n+1)/2 )^2 )/ (m1*n*(m1+n+1)/12)
  
  return(Z)
}

#####################################################################################

IC.MRL2=function(m,n,nu,p,H,per,del1,del2,del3,rho,lim,sim.num){
  set.seed(060819)
  cov.mat=matrix(rep(rho,p*p),p,p)+diag(rep((1-rho),p))
  
  RL_Array=rep(0,sim.num)
  
  for(i in 1:sim.num){
    r=matrix(0:0,nrow=m,ncol=6)
    r[,1]=1
    ref=rmvt(m,sigma=cov.mat,df=3)+r
    
    
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
      plot_stat = IPLep2(ref0,ref1,test)
    }
    if((i%%3000)==0) print(i)
    RL_Array[i]=k
  }
  obsARL<-sum(RL_Array)/sim.num # observed average run length
  sdN=sd(RL_Array) # observed sd of run length
  quantsN<-quantile(RL_Array,probs=c(0.05,0.25,0.5,0.75,0.95)) # Percentiles of run length distribution
  print(c(m,n,nu,p,H,per,rho,del1,del2,del3,obsARL,sdN,as.numeric(quantsN)))
}

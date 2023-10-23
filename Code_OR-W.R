library(lsa)
library(depth)
library(mvtnorm)
library(lcmix)

rm(list=ls())

IC.MRL2=function(m,n,p,H,per,del1,del2,del3,rho,lim,sim.num){
  set.seed(060819)
  cov.mat=matrix(rep(rho,p*p),p,p)+diag(rep((1-rho),p))
  
  RL_Array=rep(0,sim.num)
  
  for(i in 1:sim.num){
    r=matrix(0:0,nrow=m,ncol=6)
    r[,1]=1
    ref=rmvt(m,sigma=cov.mat,df=3)+r
    
    m = nrow(ref)
    
    plot_stat=k=0
    while(plot_stat<H && k<=lim){
      k=k+1
      r1=matrix(0:0,nrow=n,ncol=6)
      r1[,1]=del1
      r1[,2]=del2
      r1[,3]=r1[,4]=r1[,5]=r1[,6]=del3
      
      test=rmvt(n,sigma=cov.mat,df=3)+r1
      n = nrow(test)
      U_ref=sapply(1:m,function(i){sqrt(sum(ref[i,]^2))})
      V_test=sapply(1:n,function(i){sqrt(sum(test[i,]^2))})
      
      
      plot_stat = ((sum(rank(c(U_ref,V_test))[(m+1):(m+n)]) - n*(m+n+1)/2 )^2 )/ (m*n*(m+n+1)/12)
      
    }
    if((i%%3000)==0) print(i)
    RL_Array[i]=k
  }
  obsARL<-sum(RL_Array)/sim.num # observed average run length
  sdN=sd(RL_Array) # observed sd of run length
  quantsN<-quantile(RL_Array,probs=c(0.05,0.25,0.5,0.75,0.95)) # Percentiles of run length distribution
  print(c(m,n,rho,p,H,per,del1,del2,del3,obsARL,sdN,as.numeric(quantsN)))
}

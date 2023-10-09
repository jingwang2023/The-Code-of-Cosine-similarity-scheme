library(lsa)
library(depth)
library(mvtnorm)
library(robustbase)
library(matlib)
library(lcmix)

rm(list=ls())

OOC.MRL=function(m,n,p,H,per,del1,del2,del3,rho,lim,sim.num){
  set.seed(060819)
  cov.mat=matrix(rep(rho,p*p),p,p)+diag(rep((1-rho),p))
  RL_Array=rep(0,sim.num)
  
  
  for(u in 1:sim.num){
    
    r=matrix(0:0,nrow=m,ncol=6)
    r[,1]=1
    ref=rmvt(m,sigma=cov.mat,df=3)+r
    X=ref
    
    plot_stat=k=0
    while(plot_stat<H && k<=lim){
      
      k=k+1
      r1=matrix(0:0,nrow=n,ncol=6)
      r1[,1]=del1
      r1[,2]=del2
      r1[,3]=r1[,4]=r1[,5]=r1[,6]=del3
      test=rmvt(n,sigma=cov.mat,df=3)+r1
      Y=test
      med=colMedians(X)
      n=nrow(Y)           
      
      ### Compute the MNS charting statistic T_S^2
      p=ncol(Y)
      Vhat=array(rep(0,p*p),dim=c(p,p))
      
      for (i in 1:p){
        x[i] = sum(sign(Y[1:n,i]-med[i]))
      }
      
      for (j in 1:p){
        for (t in j:p){
          Vhat[j,t]=sum(sign(Y[1:n,j]-med[j])*sign(Y[1:n,t]-med[t]))
          Vhat[t,j]=Vhat[j,t]
          Vhat[i,i]=n
        }
      }
      
      TS2 = t(x)%*%Ginv(Vhat)%*%c(x)
      plot_stat=TS2
    }
    RL_Array[u]=k
  } 
  obsARL<-sum(RL_Array)/sim.num # observed average run length
  sdN=sd(RL_Array) # observed sd of run length
  quantsN<-quantile(RL_Array,probs=c(0.05,0.25,0.5,0.75,0.95)) # Percentiles of run length distribution
  print(c(m,n,p,H,per,del1,del2,del3,rho,obsARL,sdN,as.numeric(quantsN)))
}
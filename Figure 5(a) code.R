library(lsa)
library(depth)
library(mvtnorm)
library(MultiRNG)
library(lcmix)
library(corpcor)
library(ggplot2)

ref=read.table("D:/driver_data/train_normal.txt")
unor=read.table("D:/driver_data/test_aggressive.txt")
m=nrow(ref)
nu=round(m/5)
index = c(sample(1:m,nu))
ref0=ref[index,]
ref1=ref[-index,]

for(i in 1:27){
  
  test= unor[(30*i-29):(30*i),]
  
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
  
  Z_M = ((sum(rank(c(U_ref,V_test))[(m1+1):(m1+n)]) - n*(m1+n+1)/2 )^2 )/ (m1*n*(m1+n+1)/12)
  
  
  U1_ref=sapply(1:m1,function(i){acos(cosine(X_b_med, c(t(ref1[i,]))))})
  V1_test=sapply(1:n,function(i){acos(cosine(X_b_med, c(t(test[i,]))))})
  
  
  Z_A = ((sum(rank(c(U1_ref,V1_test))[(m1+1):(m1+n)]) - n*(m1+n+1)/2 )^2 )/ (m1*n*(m1+n+1)/12)
  
  Z=max(Z_M,Z_A)
  
  plot_stat[i]=Z  
}


x=c(1:27)
p=ggplot()+geom_line(aes(x=x,y=plot_stat),color="darkblue",size=0.5) 
p+geom_hline(aes(yintercept=9.89),color="red",size=0.5)
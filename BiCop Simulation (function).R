# cut--新息函数数据长度
# nsim--模拟次数
# aphla--模拟路径的显著性水平
# fs--Copula族范围
# S--是否进行模拟

BivariateCop=function(u1,u2,cut,nsim,fs=NA){
  
  
  ## 理论版本-----------------------------------------------------------------------
  ## Select copula and theoretical CF(concentration function)
  C=BiCopSelect(u1,u2,familyset = fs, selectioncrit = "BIC")
  
  Lth=function(z) BiCopCDF(z,z,C)/z
  Rth=function(z) (1-2*z+BiCopCDF(z,z,C))/(1-z)
  
  len=cut
  s <- seq(0, 1, length.out = len)
  
  ll_th=Vectorize(Lth)(s[1:(len/2)])
  rr_th=Vectorize(Rth)(s[((len/2)+1):len])

  lr_th=c(ll_th,rr_th)
  
  plot(x=s,y=lr_th,type='l',col='blue',xlim = c(0, 1),ylim = c(0, 1))
  
  
## 经验版本-----------------------------------------------------------------------
  U=rank(as.numeric(u1))/(nrow(u1)+1)
  V=rank(as.numeric(u2))/(nrow(u2)+1)

  # U=u1 ##
  # V=u2 ##是否要用???????????????

  Lemp=function(z) sum((U<=z)&(V<=z))/sum(U<=z)
  Remp=function(z) sum((U>=1-z)&(V>=1-z))/sum(U>=1-z)

  ll_emp=Vectorize(Lemp)(s[1:(len/2)])
  rr_emp=Vectorize(Remp)(rev(s[1:(len/2)]))

  lr_emp=c(ll_emp,rr_emp)

  lines(x=s,y=lr_emp,col='black')
  
  
  
#Simulation---------------------------------------------------------------
  
  

  MGS=matrix(NA,nsim,length(s))
  pb <- progress_bar$new(total = nsim)
  for(k in 1:nsim){
    Xs=BiCopSim(len, C)
    Us=rank(Xs[,1])/(nrow(Xs)+1)
    Vs=rank(Xs[,2])/(nrow(Xs)+1)
    Lemp=function(z) sum((Us<=z)&(Vs<=z))/sum(Us<=z)
    Remp=function(z) sum((Us>=1-z)&(Vs>=1-z))/sum(Us>=1-z)
    MGS[k,1:(length(s)/2)]=Vectorize(Lemp)(s[1:(len/2)])
    MGS[k,(length(s)/2+1):(length(s))]=Vectorize(Remp)(rev(s[1:(len/2)]))
    lines(x=s,y=MGS[k,],col="grey")
    Sys.sleep(0.01)
    pb$tick()
  }
  pb$terminate()

  
  V95=function(x) quantile(x,0.95,na.rm = T)
  Q95=apply(MGS,2,V95)
  lines(x=s,Q95,col="red",lwd=2)
  V05=function(x) quantile(x,0.05,na.rm = T)
  Q05=apply(MGS,2,V05)
  lines(x=s,Q05,col="red",lwd=2)
  
  simdf=as.data.frame(t(MGS))
  
  
  ToExcel <- cbind(s,simdf,Q95,Q05,lr_emp,lr_th)
  
  return(ToExcel)
} 
  


# a--显著性水平
# par & par2--Copula参数
# dof--自由度
# dist--残差分布
# type--Copula类型

CalculateCoVaR=function(r1,r2){

  par(mfrow=c (2,2))
  
  source("C:\\Users\\徐浩然\\Desktop\\github\\RCoVaRCopula-master\\CoVaR.R") 
  source("C:\\Users\\徐浩然\\Desktop\\github\\RCoVaRCopula-master\\CoVaRX.R") 
  source("C:\\Users\\徐浩然\\Desktop\\github\\RCoVaRCopula-master\\DynCopulaCoVaR.R") 
  source("C:\\Users\\徐浩然\\Desktop\\github\\RCoVaRCopula-master\\DynCopulaCoVaRUpper.R") 
  source("C:\\Users\\徐浩然\\Desktop\\github\\RCoVaRCopula-master\\skewtdis_inv.R") 
  require("pracma") 
  require("copula")
  
  
  ARMAGARCH=function(r,model='gjrGARCH',dist='sstd'){
    order=summary(auto.arima(r))
    AR= order[["arma"]][1]
    MA= order[["arma"]][2]
    uspec = ugarchspec(variance.model = list(model=model, garchOrder=c(1,1)),
                       mean.model=list(armaOrder=c(AR,MA), include.mean=TRUE),distribution.model=dist)
    fit = ugarchfit(data=r, spec=uspec)
    sr=rugarch::residuals(fit,standardize=T)
    uu=psstd(sr,0,1,nu = fit@fit[["coef"]][["shape"]], xi = fit@fit[["coef"]][["skew"]])
    
    return(list(u=uu,nu = fit@fit[["coef"]][["shape"]], xi = fit@fit[["coef"]][["skew"]],
                mu=fitted(fit),sig=sigma(fit)))
  }
  
  AG1=ARMAGARCH(r1)
  AG2=ARMAGARCH(r2)
  
  CC=BiCopSelect(AG1[["u"]],AG2[["u"]],familyset = 2, selectioncrit = "BIC")
  print(CC)
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # VaR & CVaR
  Q1=CoVaRX(0.05,0.05,par=CC[["par"]],par2=CC[["par2"]],dof=AG1[["nu"]],gamma=AG1[["xi"]],cond.mean=AG1[["mu"]], cond.sigma=AG1[["sig"]])
  BM1=CoVaRX(0.05,0.5,par=CC[["par"]],par2=CC[["par2"]],dof=AG1[["nu"]],gamma=AG1[["xi"]],cond.mean=AG1[["mu"]], cond.sigma=AG1[["sig"]])
  VaRD1=Q1$VaR
  VaRU1=Q1$VaR.up
  CoVaRD1=Q1$CoVaR
  CoVaRU1=Q1$CoVaR.up
  CBMD1=BM1$CoVaR
  CBMU1=BM1$CoVaR.up

  # delta CoVaR
  deltaCDw1=CoVaRD1-CBMD1
  deltaCUp1=CoVaRU1-CBMU1
  deltaCDw1_pct=(CoVaRD1-CBMD1)/CBMD1
  deltaCUp1_pct=(CoVaRU1-CBMU1)/CBMU1
  
  # plot
  plot(VaRD1,type='l',col='blue',ylim = c(min(CoVaRD1), max(CoVaRU1)))
  lines(VaRU1,col='blue')
  lines(CoVaRD1,col='black')
  lines(CoVaRU1,col='black')
  lines(CBMD1,col='red')
  lines(CBMU1,col='red')
  lines(deltaCDw1,col='green')
  lines(deltaCUp1,col='green')
  
  plot(deltaCDw1_pct,type='l',col='red',ylim = c(min(deltaCDw1_pct,deltaCUp1_pct), max(deltaCDw1_pct,deltaCUp1_pct)))
  lines(deltaCUp1_pct,col='blue')
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  
  # VaR & CVaR
  Q2=CoVaRX(0.05,0.05,par=CC[["par"]],par2=CC[["par2"]],dof=AG2[["nu"]],gamma=AG2[["xi"]],cond.mean=AG2[["mu"]], cond.sigma=AG2[["sig"]])
  BM2=CoVaRX(0.05,0.5,par=CC[["par"]],par2=CC[["par2"]],dof=AG2[["nu"]],gamma=AG2[["xi"]],cond.mean=AG2[["mu"]], cond.sigma=AG2[["sig"]])
  VaRD2=Q2$VaR
  VaRU2=Q2$VaR.up
  CoVaRD2=Q2$CoVaR
  CoVaRU2=Q2$CoVaR.up
  CBMD2=BM2$CoVaR
  CBMU2=BM2$CoVaR.up
  
  # delta CoVaR
  deltaCDw2=CoVaRD2-CBMD2
  deltaCUp2=CoVaRU2-CBMU2
  deltaCDw2_pct=(CoVaRD2-CBMD2)/CBMD2
  deltaCUp2_pct=(CoVaRU2-CBMU2)/CBMU2
  
  # plot
  plot(VaRD2,type='l',col='blue',ylim = c(min(CoVaRD2), max(CoVaRU2)))
  lines(VaRU2,col='blue')
  lines(CoVaRD2,col='black')
  lines(CoVaRU2,col='black')
  lines(CBMD2,col='red')
  lines(CBMU2,col='red')
  lines(deltaCDw2,col='green')
  lines(deltaCUp2,col='green')
  
  plot(deltaCDw2_pct,type='l',col='red',ylim = c(min(deltaCDw2_pct,deltaCUp2_pct), max(deltaCDw2_pct,deltaCUp2_pct)))
  lines(deltaCUp2_pct,col='blue')
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  print("Is significant? 2->1 D&U")
  print(ks.test(CoVaRD1, VaRD1))
  print(ks.test(CoVaRU1, VaRU1))
  
  print("Is significant? 1->2 D&U")
  print(ks.test(CoVaRD2, VaRD2))
  print(ks.test(CoVaRU2, VaRU2))
  
  print("Asymmtric? 1 & 2")
  print(ks.test(deltaCDw1_pct,   deltaCUp1_pct))
  print(ks.test(deltaCDw2_pct,   deltaCUp2_pct))
  
  
  
  par(mfrow=c(1, 1)) 
  
  TOEXL=as.data.frame(cbind(VaRD1,VaRU1,CoVaRD1,CoVaRU1,CBMD1,CBMU1,deltaCDw1,deltaCUp1,deltaCDw1_pct,deltaCUp1_pct,
                            VaRD2,VaRU2,CoVaRD2,CoVaRU2,CBMD2,CBMU2,deltaCDw2,deltaCUp2,deltaCDw2_pct,deltaCUp2_pct))
  TOEXL=cbind(dates_1,TOEXL)
  return (TOEXL)
}


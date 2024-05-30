
## r--单个对数收益率
## alpha--计算VaR的分位数水平
## newslim--新息的x轴范围
## model--GARCH模型
## dist--GARCH模型的数据分布
##记得修改ARMA order部分
##
## Note:model Valid models (currently implemented) are 
## “sGARCH”, “fGARCH”, “eGARCH”, “gjrGARCH”, “apARCH” and “iGARCH” and “csGARCH”.
## submodel If the model is “fGARCH”, valid submodels are 
##“GARCH”, “TGARCH”, “AVGARCH”, “NGARCH”, “NAGARCH”, “APARCH”,“GJRGARCH” and “ALLGARCH”
##
## Valid choices are “norm” for the normal distibution, 
## “snorm” for the skew-normal distribution, 
## “std” for the student-t, “sstd” for the skew-student, 
## “ged” for the generalized error distribution, 
## “sged” for the skew-generalized error distribution, 
## “nig” for the normal inverse gaussian distribution, “ghyp”

ARMA_GRACH=function(r,sn,alpha,newslim,model,dist,lag){

  
order=summary(auto.arima(r,d=0,seasonal = F,ic='bic',stepwise = F))
AR= order[["arma"]][1]#!!!!!!!!!!!!!!!!
MA= order[["arma"]][2]#!!!!!!!!!!!!!!!!


# 1: 3,3
uspec = ugarchspec(variance.model = list(model=model, garchOrder=c(1,1)),
                    mean.model=list(armaOrder=c(AR,MA), include.mean=F),distribution.model=dist)##分布的选择

fit = ugarchfit(data=r, spec=uspec)

sink(paste("C:\\Users\\徐浩然\\Desktop\\P4 Data\\O-",sn,"_armagarch.txt",sep = ""),split=T)
print(fit)
sink()

mu=fitted(fit)   ##到底是fitted or uncmean??????????????????????????????????????????
sig=sigma(fit)

sre=rugarch::residuals(fit,standardize=T)

print(ks.test(sre,'psstd',0,1,nu = fit@fit[["coef"]][["shape"]], xi = fit@fit[["coef"]][["skew"]]))#----------------------注意

u=psstd(sre,0,1,nu = fit@fit[["coef"]][["shape"]], xi = fit@fit[["coef"]][["skew"]])


#News impact curve
ni=newsimpact(fit, seq(-newslim, newslim, length.out = length(r)))
newsx=ni$zx
newsy=ni$zy








# Calculate VaR and plot data with VaR
VaR = function(x, a)
{
  vmodel  = x@model$modeldesc$vmodel
  T = x@model$modeldata$T
  insample = 1:T
  xseries = x@model$modeldata$data[insample]
  xdates  = x@model$modeldata$index[insample]
  xsigma 	= x@fit$sigma
  distribution = x@model$modeldesc$distribution
  xcmu = fitted(x)
  idx = x@model$pidx
  pars  = x@fit$ipars[,1]
  skew  = pars[idx["skew",1]]
  shape = pars[idx["shape",1]]
  if(distribution == "ghst") ghlambda = -shape/2 else ghlambda = pars[idx["ghlambda",1]]
  z1 	= a
  z2 	= 1-a
  
  qDw 	= fitted(x) + sigma(x)* qdist(distribution, z1, 0, 1, lambda = ghlambda, skew, shape)
  qUp	= fitted(x) + sigma(x)* qdist(distribution, z2, 0, 1, lambda = ghlambda, skew, shape)
  
# Plot-------------------------------------------------------------------------------------
  plot(xdates, xseries, type = "l", col = "black", ylab = "Returns", xlab="Time",ylim = c(-25, 25),
       main = paste("Series with with",a*100,"% VaR Limits"), cex.main = 0.8)
  lines(xdates, qDw, col = "red",lwd=2,)
  lines(xdates, qUp, col = "blue",lwd=2,)
  mtext(paste("GARCH model :", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
  if(vmodel == "fGARCH"){
    mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
  }
  abline(h = 0, col = "grey", lty = 3)
  grid()
# Plot-------------------------------------------------------------------------------------
}

VaR=VaR(fit,alpha)
VaR_Dw=VaR$qDw
VaR_Up=VaR$qUp

return (cbind(dates_1,as.data.frame(cbind(mu,sig,u,r,sre,VaR_Dw,VaR_Up,newsx,newsy))))
}








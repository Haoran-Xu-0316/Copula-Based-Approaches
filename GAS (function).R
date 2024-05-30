# u--残差的数据框
# dist--残差服从的分布

GAS=function(u,dist,location=F,scale=F,correlation=T,shape=T){

GASSpec = MultiGASSpec(Dist = dist, ScalingType = "Identity",
             GASPar = list(location=location,scale = scale, correlation = correlation,shape=shape))
Fit = MultiGASFit(GASSpec, u)
summary(Fit)
plot(Fit,which=1)# can be deleted
par(mfrow = c(1, 1))
getPar=as.data.frame(getFilteredParameters(Fit))
# rho=getPar$rho12

return (cbind(dates,getPar))
}







tailtau=function(u1,u2,fs=NA,par,par2){
C=BiCopSelect(u1,u2,familyset = 2)
print(C)
# Time-varying tau and tail dependence
taildep = BiCopPar2TailDep(family = 2#C[["family"]],
                           ,par=par,par2=cop$par2)
uptail=taildep$upper
lowtail=taildep$lower

tau = BiCopPar2Tau(family = 2, par = par,par2=cop$par2)

plot(tau,type='l')
return (cbind(dates,tau,uptail,lowtail))
}


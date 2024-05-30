library(TSP)
library(xts)
library(GAS)
library(aTSA)
library(vars)
library(stats)
library(tidyr)
library(FinTS)
library(dplyr)
library(vines)
library(readxl)
library(pracma)
library(fGarch)
library(copula)
library(tseries)
library(ggplot2)
library(rugarch)
library(rmgarch)
library(tseries)
library(moments)
library(forecast)
library(openxlsx)
library(progress)
library(VineCopula)
library(rvinecopulib)
library(PerformanceAnalytics)
library(CDVineCopulaConditional)


# path--数据文件路径
# sheet--excel文件的工作表

read_data=function(path,sheet)

{
  
data = read_excel(path, sheet = sheet)
label=colnames(data)
label=label[-1]
dates <- as.Date(data$Date, format = "%Y%m%d")
data <- xts(data[, label], order.by = dates)

rt = 100*apply(log(data),2,diff)
write.csv(rt, file="C:\\Users\\徐浩然\\Desktop\\P4 Data\\rt.csv")
return (list(rt,dates))

}

# 
stat=function(data,R=F){
  if (R==T){  
    
    Series     = colnames(data)
    N          = length(Series)
    stats      = matrix(0,6,N,dimnames=list(c('L-B', 'p-value' ,'L-B2','p-value','ARCH','p-value'),Series))
    
    for (i in 1:N){
      
      stats[1,i] = Box.test(data[,Series[i]],10,"Ljung-Box")$statistic
      stats[2,i] = Box.test((data[,Series[i]])^2,10, "Ljung-Box")$p.value
      stats[3,i] = Box.test((data[,Series[i]])^2,10,"Ljung-Box")$statistic
      stats[4,i] = Box.test(data[,Series[i]],10, "Ljung-Box")$p.value
      stats[5,i] = ArchTest (data[,Series[i]],10)$statistic
      stats[6,i] = ArchTest (data[,Series[i]],10)$p.value
    }
    print(t(stats), digits=6)
    write.csv(stats, file="C:\\Users\\徐浩然\\Desktop\\P4 Data\\Residual.csv")}
  
  
else{
  Series     = colnames(data)
  N          = length(Series)
  stats      = matrix(0,14,N,dimnames=list(c('Mean','Max.','Min.','Std. Dev.','Skewness','Kurtosis','J-B', 'p-value' ,'ADF','p-value','L-B','p-value','ARCH','p-value'),Series))

  for (i in 1:N){
    lg=20
    stats[1,i]  = mean(data[,Series[i]])
    stats[2,i]  = max(data[,Series[i]])
    stats[3,i]  = min(data[,Series[i]])
    stats[4,i]  = sd(data[,Series[i]])
    stats[5,i]  = skewness(data[,Series[i]])
    stats[6,i]  = kurtosis(data[,Series[i]])
    stats[7,i]  = jarque.bera.test(data[,Series[i]])$statistic
    stats[8,i]  = jarque.bera.test(data[,Series[i]])$p.value
    stats[9,i] = tseries::adf.test(data[,Series[i]],k = 10)$statistic
    stats[10,i] = tseries::adf.test(data[,Series[i]],k = 10)$p.value
    stats[11,i] = Box.test(data[,Series[i]],lg,"Ljung-Box")$statistic
    stats[12,i] = Box.test(data[,Series[i]],lg, "Ljung-Box")$p.value
    stats[13,i] = ArchTest (data[,Series[i]],lg)$statistic
    stats[14,i] = ArchTest (data[,Series[i]],lg)$p.value
  }
  print(t(stats), digits=6)
  write.csv(stats, file="C:\\Users\\徐浩然\\Desktop\\P4 Data\\Descriptive.csv")
}
}


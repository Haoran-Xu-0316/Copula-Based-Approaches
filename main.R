options(warn=-1)

source('C:\\Users\\徐浩然\\Desktop\\R file\\Paper 4\\Function\\GAS (function).R')
source('C:\\Users\\徐浩然\\Desktop\\R file\\Paper 4\\Function\\CoVaR (function).R')
source('C:\\Users\\徐浩然\\Desktop\\R file\\Paper 4\\Function\\ReadData (function).R')
source('C:\\Users\\徐浩然\\Desktop\\R file\\Paper 4\\Function\\ARMA-GARCH (function).R')
source('C:\\Users\\徐浩然\\Desktop\\R file\\Paper 4\\Function\\VineCopula (function).R')
source('C:\\Users\\徐浩然\\Desktop\\R file\\Paper 4\\Function\\BiCop Simulation (function).R')



# Modeling processes
rd=read_data("C:\\Users\\徐浩然\\Desktop\\P4 Data\\Sectordata.xlsx","Sheet1")
dates_1=rd[[2]][-1]
dates=rd[[2]]
rt=as.data.frame(rd[1])
sheets=colnames(rt)
# x=cor(rt)
stat(rt)
## 11 NESG 12 ESG
W=12
# ARMA-GARCH modeling
aa=0.05
ns=0.3
dt="sstd"
m="gjrGARCH"
step1_1 = ARMA_GRACH(rt[,1],sn=sheets[1],a=aa,newslim=ns,dist=dt,model=m,lag=c(3,3))
step1_2 = ARMA_GRACH(rt[,2],sn=sheets[2],a=aa,newslim=ns,dist=dt,model=m,lag=c(3,3))
step1_3 = ARMA_GRACH(rt[,3],sn=sheets[3],a=aa,newslim=ns,dist=dt,model=m,lag=c(3,3))
step1_4 = ARMA_GRACH(rt[,4],sn=sheets[4],a=aa,newslim=ns,dist=dt,model=m,lag=c(2,3))
step1_5 = ARMA_GRACH(rt[,5],sn=sheets[5],a=aa,newslim=ns,dist=dt,model=m,lag=c(2,3))
step1_6 = ARMA_GRACH(rt[,6],sn=sheets[6],a=aa,newslim=ns,dist=dt,model=m,lag=c(2,2))
step1_7 = ARMA_GRACH(rt[,7],sn=sheets[7],a=aa,newslim=ns,dist=dt,model=m,lag=c(2,2))
step1_8 = ARMA_GRACH(rt[,8],sn=sheets[8],a=aa,newslim=ns,dist=dt,model=m,lag=c(3,2))
step1_9 = ARMA_GRACH(rt[,9],sn=sheets[9],a=aa,newslim=ns,dist=dt,model=m,lag=c(2,3))
step1_10 = ARMA_GRACH(rt[,10],sn=sheets[10],a=aa,newslim=ns,dist=dt,model=m,lag=c(2,2))
step1_11 = ARMA_GRACH(rt[,11],sn=sheets[11],a=aa,newslim=ns,dist=dt,model=m,lag=c(2,2))
step1_12 = ARMA_GRACH(rt[,12],sn=sheets[12],a=aa,newslim=ns,dist=dt,model=m,lag=c(2,2))





# Extracting 
u=data.frame(step1_1$u,step1_2$u,step1_3$u,step1_4$u,step1_5$u,step1_6$u,
             step1_7$u,step1_8$u,step1_9$u,step1_10$u,step1_11$u,step1_12$u)
sre=data.frame(step1_1$sre,step1_2$sre,step1_3$sre,step1_4$sre,step1_5$sre,step1_6$sre,
               step1_7$sre,step1_8$sre,step1_9$sre,step1_10$sre,step1_11$sre,step1_12$sre)
u <- xts(u, order.by = dates_1)
sre <- xts(sre, order.by = dates_1)
colnames(u)=sheets
colnames(sre)=sheets

stat(sre,R=T)


for (i in 1:10){
  cop=BiCopSelect(u[,i],u[,12],familyset=2)
  print(cop)
  print(cop$taildep$lower)
  # x=BiCopGofTest(u[,i],u[,12],family=2,method = "white",B = 50)
  # print(x)

}
# GAS modeling 
lo=T
sl=F
pho=T
sp=F
d='mvt'## mvnorm
step3_1=GAS(u[,c(1,W)],dist=d,location=lo,scale=sl,correlation=pho,shape=sp)
step3_2=GAS(u[,c(2,W)],dist=d,location=lo,scale=sl,correlation=pho,shape=sp)
step3_3=GAS(u[,c(3,W)],dist=d,location=lo,scale=sl,correlation=pho,shape=sp)
step3_4=GAS(u[,c(4,W)],dist=d,location=lo,scale=sl,correlation=pho,shape=sp)
step3_5=GAS(u[,c(5,W)],dist=d,location=lo,scale=sl,correlation=pho,shape=sp)
step3_6=GAS(u[,c(6,W)],dist=d,location=lo,scale=sl,correlation=pho,shape=sp)
step3_7=GAS(u[,c(7,W)],dist=d,location=lo,scale=sl,correlation=pho,shape=sp)
step3_8=GAS(u[,c(8,W)],dist=d,location=lo,scale=sl,correlation=pho,shape=sp)
step3_9=GAS(u[,c(9,W)],dist=d,location=lo,scale=sl,correlation=pho,shape=sp)
step3_10=GAS(u[,c(10,W)],dist=d,location=lo,scale=sl,correlation=pho,shape=sp)






#  Extract tau & tail dependence
tau_tail_1=tailtau(u[,1],u[,W],fs=NA,par=step3_1$rho12,par2=step3_1$nu)
tau_tail_2=tailtau(u[,2],u[,W],fs=NA,par=step3_2$rho12,par2=step3_2$nu)
tau_tail_3=tailtau(u[,3],u[,W],fs=NA,par=step3_3$rho12,par2=step3_3$nu)
tau_tail_4=tailtau(u[,4],u[,W],fs=NA,par=step3_4$rho12,par2=step3_4$nu)
tau_tail_5=tailtau(u[,5],u[,W],fs=NA,par=step3_5$rho12,par2=step3_5$nu)
tau_tail_6=tailtau(u[,6],u[,W],fs=NA,par=step3_6$rho12,par2=step3_6$nu)
tau_tail_7=tailtau(u[,7],u[,W],fs=NA,par=step3_7$rho12,par2=step3_7$nu)
tau_tail_8=tailtau(u[,8],u[,W],fs=NA,par=step3_8$rho12,par2=step3_8$nu)
tau_tail_9=tailtau(u[,9],u[,W],fs=NA,par=step3_9$rho12,par2=step3_9$nu)
tau_tail_10=tailtau(u[,10],u[,W],fs=NA,par=step3_10$rho12,par2=step3_10$nu)



# Calculating CoVaR
step4_1=CalculateCoVaR(rt[,1],rt[,W])
step4_2=CalculateCoVaR(rt[,2],rt[,W])
step4_3=CalculateCoVaR(rt[,3],rt[,W])
step4_4=CalculateCoVaR(rt[,4],rt[,W])
step4_5=CalculateCoVaR(rt[,5],rt[,W])
step4_6=CalculateCoVaR(rt[,6],rt[,W])
step4_7=CalculateCoVaR(rt[,7],rt[,W])
step4_8=CalculateCoVaR(rt[,8],rt[,W])
step4_9=CalculateCoVaR(rt[,9],rt[,W])
step4_10=CalculateCoVaR(rt[,10],rt[,W])




# Vine-Copula modeling
step2_R=vineCopula(u[,1:10],"R",fs=NA,selectc="AIC")
step2_C=vineCopula(u[,1:10],"C",fs=NA,selectc="AIC")
step2_D=vineCopula(u[,1:10],"D",fs=NA,selectc="AIC")



RVineVuongTest(u[,1:10], step2_R, step2_C)
RVineVuongTest(u[,1:10], step2_R, step2_D)
RVineVuongTest(u[,1:10], step2_C, step2_D)


RVineClarkeTest(u[,1:10], step2_R, step2_C)
RVineClarkeTest(u[,1:10], step2_R, step2_D)
RVineClarkeTest(u[,1:10], step2_C, step2_D)


cut=500
nsim=150
fa=2
# Simulation of concentration function
step5_1=BivariateCop(u[,1],u[,W],cut=cut,nsim=nsim,fs=fa)
step5_2=BivariateCop(u[,2],u[,W],cut=cut,nsim=nsim,fs=fa)
step5_3=BivariateCop(u[,3],u[,W],cut=cut,nsim=nsim,fs=fa)
step5_4=BivariateCop(u[,4],u[,W],cut=cut,nsim=nsim,fs=fa)
step5_5=BivariateCop(u[,5],u[,W],cut=cut,nsim=nsim,fs=fa)
step5_6=BivariateCop(u[,6],u[,W],cut=cut,nsim=nsim,fs=fa)
step5_7=BivariateCop(u[,7],u[,W],cut=cut,nsim=nsim,fs=fa)
step5_8=BivariateCop(u[,8],u[,W],cut=cut,nsim=nsim,fs=fa)
step5_9=BivariateCop(u[,9],u[,W],cut=cut,nsim=nsim,fs=fa)
step5_10=BivariateCop(u[,10],u[,W],cut=cut,nsim=nsim,fs=fa)











# Save data to excel
if (W==12){
wb0 <- createWorkbook()
addWorksheet(wb0, 'U')
writeData(wb0, 'U', cbind(dates_1,as.data.frame(u)))
addWorksheet(wb0, 'rt')
writeData(wb0, 'rt', cbind(dates_1,rt))
addWorksheet(wb0, 'sre')
writeData(wb0, 'sre', cbind(dates_1,as.data.frame(sre)))
saveWorkbook(wb0, "C:\\Users\\徐浩然\\Desktop\\P4 ESG\\ReadData.xlsx", overwrite = TRUE)



wb1 <- createWorkbook()
addWorksheet(wb1, sheets[1])
writeData(wb1, sheets[1], step1_1)
addWorksheet(wb1, sheets[2])
writeData(wb1, sheets[2], step1_2)
addWorksheet(wb1, sheets[3])
writeData(wb1, sheets[3], step1_3)
addWorksheet(wb1, sheets[4])
writeData(wb1, sheets[4], step1_4)
addWorksheet(wb1, sheets[5])
writeData(wb1, sheets[5], step1_5)
addWorksheet(wb1, sheets[6])
writeData(wb1, sheets[6], step1_6)
addWorksheet(wb1, sheets[7])
writeData(wb1, sheets[7], step1_7)
addWorksheet(wb1, sheets[8])
writeData(wb1, sheets[8], step1_8)
addWorksheet(wb1, sheets[9])
writeData(wb1, sheets[9], step1_9)
addWorksheet(wb1, sheets[10])
writeData(wb1, sheets[10], step1_10)
saveWorkbook(wb1, "C:\\Users\\徐浩然\\Desktop\\P4 ESG\\ARMA_GARCH.xlsx", overwrite = TRUE)


wb3 <- createWorkbook()
addWorksheet(wb3, sheets[1])
writeData(wb3, sheets[1], step3_1)
addWorksheet(wb3, sheets[2])
writeData(wb3, sheets[2], step3_2)
addWorksheet(wb3, sheets[3])
writeData(wb3, sheets[3], step3_3)
addWorksheet(wb3, sheets[4])
writeData(wb3, sheets[4], step3_4)
addWorksheet(wb3, sheets[5])
writeData(wb3, sheets[5], step3_5)
addWorksheet(wb3, sheets[6])
writeData(wb3, sheets[6], step3_6)
addWorksheet(wb3, sheets[7])
writeData(wb3, sheets[7], step3_7)
addWorksheet(wb3, sheets[8])
writeData(wb3, sheets[8], step3_8)
addWorksheet(wb3, sheets[9])
writeData(wb3, sheets[9], step3_9)
addWorksheet(wb3, sheets[10])
writeData(wb3, sheets[10], step3_10)
saveWorkbook(wb3, "C:\\Users\\徐浩然\\Desktop\\P4 ESG\\GAS.xlsx", overwrite = TRUE)



wb33 <- createWorkbook()
addWorksheet(wb33, sheets[1])
writeData(wb33, sheets[1], tau_tail_1)
addWorksheet(wb33, sheets[2])
writeData(wb33, sheets[2], tau_tail_2)
addWorksheet(wb33, sheets[3])
writeData(wb33, sheets[3], tau_tail_3)
addWorksheet(wb33, sheets[4])
writeData(wb33, sheets[4], tau_tail_4)
addWorksheet(wb33, sheets[5])
writeData(wb33, sheets[5], tau_tail_5)
addWorksheet(wb33, sheets[6])
writeData(wb33, sheets[6], tau_tail_6)
addWorksheet(wb33, sheets[7])
writeData(wb33, sheets[7], tau_tail_7)
addWorksheet(wb33, sheets[8])
writeData(wb33, sheets[8], tau_tail_8)
addWorksheet(wb33, sheets[9])
writeData(wb33, sheets[9], tau_tail_9)
addWorksheet(wb33, sheets[10])
writeData(wb33, sheets[10], tau_tail_10)
saveWorkbook(wb33, "C:\\Users\\徐浩然\\Desktop\\P4 ESG\\Tau_Taildep.xlsx", overwrite = TRUE)


wb4 <- createWorkbook()
addWorksheet(wb4, sheets[1])
writeData(wb4, sheets[1], step4_1)
addWorksheet(wb4, sheets[2])
writeData(wb4, sheets[2], step4_2)
addWorksheet(wb4, sheets[3])
writeData(wb4, sheets[3], step4_3)
addWorksheet(wb4, sheets[4])
writeData(wb4, sheets[4], step4_4)
addWorksheet(wb4, sheets[5])
writeData(wb4, sheets[5], step4_5)
addWorksheet(wb4, sheets[6])
writeData(wb4, sheets[6], step4_6)
addWorksheet(wb4, sheets[7])
writeData(wb4, sheets[7], step4_7)
addWorksheet(wb4, sheets[8])
writeData(wb4, sheets[8], step4_8)
addWorksheet(wb4, sheets[9])
writeData(wb4, sheets[9], step4_9)
addWorksheet(wb4, sheets[10])
writeData(wb4, sheets[10], step4_10)
saveWorkbook(wb4, "C:\\Users\\徐浩然\\Desktop\\P4 ESG\\CoVaR.xlsx", overwrite = TRUE)


wb5 <- createWorkbook()
addWorksheet(wb5, sheets[1])
writeData(wb5, sheets[1], step5_1)
addWorksheet(wb5, sheets[2])
writeData(wb5, sheets[2], step5_2)
addWorksheet(wb5, sheets[3])
writeData(wb5, sheets[3], step5_3)
addWorksheet(wb5, sheets[4])
writeData(wb5, sheets[4], step5_4)
addWorksheet(wb5, sheets[5])
writeData(wb5, sheets[5], step5_5)
addWorksheet(wb5, sheets[6])
writeData(wb5, sheets[6], step5_6)
addWorksheet(wb5, sheets[7])
writeData(wb5, sheets[7], step5_7)
addWorksheet(wb5, sheets[8])
writeData(wb5, sheets[8], step5_8)
addWorksheet(wb5, sheets[9])
writeData(wb5, sheets[9], step5_9)
addWorksheet(wb5, sheets[10])
writeData(wb5, sheets[10], step5_10)
saveWorkbook(wb5, "C:\\Users\\徐浩然\\Desktop\\P4 ESG\\Concentration Function.xlsx", overwrite = TRUE)

}


if (W==11){
  wb0 <- createWorkbook()
  addWorksheet(wb0, 'U')
  writeData(wb0, 'U', cbind(dates_1,as.data.frame(u)))
  addWorksheet(wb0, 'rt')
  writeData(wb0, 'rt', cbind(dates_1,rt))
  addWorksheet(wb0, 'sre')
  writeData(wb0, 'sre', cbind(dates_1,as.data.frame(sre)))
  saveWorkbook(wb0, "C:\\Users\\徐浩然\\Desktop\\P4 NOESG\\ReadData.xlsx", overwrite = TRUE)
  
  
  
  wb1 <- createWorkbook()
  addWorksheet(wb1, sheets[1])
  writeData(wb1, sheets[1], step1_1)
  addWorksheet(wb1, sheets[2])
  writeData(wb1, sheets[2], step1_2)
  addWorksheet(wb1, sheets[3])
  writeData(wb1, sheets[3], step1_3)
  addWorksheet(wb1, sheets[4])
  writeData(wb1, sheets[4], step1_4)
  addWorksheet(wb1, sheets[5])
  writeData(wb1, sheets[5], step1_5)
  addWorksheet(wb1, sheets[6])
  writeData(wb1, sheets[6], step1_6)
  addWorksheet(wb1, sheets[7])
  writeData(wb1, sheets[7], step1_7)
  addWorksheet(wb1, sheets[8])
  writeData(wb1, sheets[8], step1_8)
  addWorksheet(wb1, sheets[9])
  writeData(wb1, sheets[9], step1_9)
  addWorksheet(wb1, sheets[10])
  writeData(wb1, sheets[10], step1_10)
  saveWorkbook(wb1, "C:\\Users\\徐浩然\\Desktop\\P4 NOESG\\ARMA_GARCH.xlsx", overwrite = TRUE)
  
  
  wb3 <- createWorkbook()
  addWorksheet(wb3, sheets[1])
  writeData(wb3, sheets[1], step3_1)
  addWorksheet(wb3, sheets[2])
  writeData(wb3, sheets[2], step3_2)
  addWorksheet(wb3, sheets[3])
  writeData(wb3, sheets[3], step3_3)
  addWorksheet(wb3, sheets[4])
  writeData(wb3, sheets[4], step3_4)
  addWorksheet(wb3, sheets[5])
  writeData(wb3, sheets[5], step3_5)
  addWorksheet(wb3, sheets[6])
  writeData(wb3, sheets[6], step3_6)
  addWorksheet(wb3, sheets[7])
  writeData(wb3, sheets[7], step3_7)
  addWorksheet(wb3, sheets[8])
  writeData(wb3, sheets[8], step3_8)
  addWorksheet(wb3, sheets[9])
  writeData(wb3, sheets[9], step3_9)
  addWorksheet(wb3, sheets[10])
  writeData(wb3, sheets[10], step3_10)
  saveWorkbook(wb3, "C:\\Users\\徐浩然\\Desktop\\P4 NOESG\\GAS.xlsx", overwrite = TRUE)
  
  
  
  wb33 <- createWorkbook()
  addWorksheet(wb33, sheets[1])
  writeData(wb33, sheets[1], tau_tail_1)
  addWorksheet(wb33, sheets[2])
  writeData(wb33, sheets[2], tau_tail_2)
  addWorksheet(wb33, sheets[3])
  writeData(wb33, sheets[3], tau_tail_3)
  addWorksheet(wb33, sheets[4])
  writeData(wb33, sheets[4], tau_tail_4)
  addWorksheet(wb33, sheets[5])
  writeData(wb33, sheets[5], tau_tail_5)
  addWorksheet(wb33, sheets[6])
  writeData(wb33, sheets[6], tau_tail_6)
  addWorksheet(wb33, sheets[7])
  writeData(wb33, sheets[7], tau_tail_7)
  addWorksheet(wb33, sheets[8])
  writeData(wb33, sheets[8], tau_tail_8)
  addWorksheet(wb33, sheets[9])
  writeData(wb33, sheets[9], tau_tail_9)
  addWorksheet(wb33, sheets[10])
  writeData(wb33, sheets[10], tau_tail_10)
  saveWorkbook(wb33, "C:\\Users\\徐浩然\\Desktop\\P4 NOESG\\Tau_Taildep.xlsx", overwrite = TRUE)
  
  
  wb4 <- createWorkbook()
  addWorksheet(wb4, sheets[1])
  writeData(wb4, sheets[1], step4_1)
  addWorksheet(wb4, sheets[2])
  writeData(wb4, sheets[2], step4_2)
  addWorksheet(wb4, sheets[3])
  writeData(wb4, sheets[3], step4_3)
  addWorksheet(wb4, sheets[4])
  writeData(wb4, sheets[4], step4_4)
  addWorksheet(wb4, sheets[5])
  writeData(wb4, sheets[5], step4_5)
  addWorksheet(wb4, sheets[6])
  writeData(wb4, sheets[6], step4_6)
  addWorksheet(wb4, sheets[7])
  writeData(wb4, sheets[7], step4_7)
  addWorksheet(wb4, sheets[8])
  writeData(wb4, sheets[8], step4_8)
  addWorksheet(wb4, sheets[9])
  writeData(wb4, sheets[9], step4_9)
  addWorksheet(wb4, sheets[10])
  writeData(wb4, sheets[10], step4_10)
  saveWorkbook(wb4, "C:\\Users\\徐浩然\\Desktop\\P4 NOESG\\CoVaR.xlsx", overwrite = TRUE)
  
  
  wb5 <- createWorkbook()
  addWorksheet(wb5, sheets[1])
  writeData(wb5, sheets[1], step5_1)
  addWorksheet(wb5, sheets[2])
  writeData(wb5, sheets[2], step5_2)
  addWorksheet(wb5, sheets[3])
  writeData(wb5, sheets[3], step5_3)
  addWorksheet(wb5, sheets[4])
  writeData(wb5, sheets[4], step5_4)
  addWorksheet(wb5, sheets[5])
  writeData(wb5, sheets[5], step5_5)
  addWorksheet(wb5, sheets[6])
  writeData(wb5, sheets[6], step5_6)
  addWorksheet(wb5, sheets[7])
  writeData(wb5, sheets[7], step5_7)
  addWorksheet(wb5, sheets[8])
  writeData(wb5, sheets[8], step5_8)
  addWorksheet(wb5, sheets[9])
  writeData(wb5, sheets[9], step5_9)
  addWorksheet(wb5, sheets[10])
  writeData(wb5, sheets[10], step5_10)
  saveWorkbook(wb5, "C:\\Users\\徐浩然\\Desktop\\P4 NOESG\\Concentration Function.xlsx", overwrite = TRUE)
  
}


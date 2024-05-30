# uda--残差的数据框
# type--vine的类型(Rvine,Dvine,Cvine)
# fs--Copula族范围
# selectc--模型选择标准(AIC,BIC...)
vineCopula=function(uda,type,fs=NA,selectc){
  
  if (type=="D"){
    d <- dim(uda)[2]
    M <- 1 - abs(TauMatrix(uda))
    hamilton <- insert_dummy(TSP(M), label = "cut")
    sol <- solve_TSP(hamilton, method = "repetitive_nn")
    order <- cut_tour(sol, "cut")
    DVM <- D2RVine(order, family = rep(0,d*(d-1)/2), par = rep(0, d*(d-1)/2))
    RV=RVineCopSelect(uda,familyset = fs, Matrix=DVM$Matrix, selectioncrit = selectc,rotations=T)
    # RV=RVineStructureSelect(uda, familyset = fs, type = 0, selectioncrit = selectc,rotations=T)
  }
  
  if(type=="C"){RV=RVineStructureSelect(uda, familyset = fs, type = 1, selectioncrit = selectc,rotations=T)}
  
  if(type=="R"){RV=RVineStructureSelect(uda, familyset = fs, type = 0, selectioncrit = selectc,rotations=T)}
  
  
  sink(paste("C:\\Users\\徐浩然\\Desktop\\P4 Data\\",type,"vine.txt",sep = ""),split=T)
  summary(RV)
  sink()
  plot(RV, tree = "ALL", type = 1,edge.labels='family-tau')

# RVineTreePlot(
#   RV,
#   tree = "ALL",
#   type = 0,
#   edge.labels = NULL,
#   legend.pos = "bottomleft",
#   interactive = T,
#   ...
# )
  return(RV)
}


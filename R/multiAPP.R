#' @export
multiAPP=function(HA,PS,mc){
  # multivariate approach
  GS=parallel::mclapply(1:length(HA$Y),function(c){
    return(GCN(HA$S,as.matrix(HA$M[[c]]),HA$Y[[c]],as.data.frame(HA$X),PS[[c]]))
  },mc.cores = mc,mc.preschedule = T,mc.cleanup = T)
  # GS=foreach::foreach(c=1:length(HA$Y)) %dopar% { GCN(HA$S,HA$M[[c]],HA$Y[[c]],HA$X,PS[[c]]) }
  names(GS)=names(HA$Y)
  Q=lapply(GS,function(gs){ return(data.table::rbindlist(gs,idcol = "app")) })
  QQ=data.table::rbindlist(Q,idcol = "ensg")
  GP=data.table::rbindlist(lapply(split(QQ,QQ$app),ppp,mode = "p"))
  GG=cbind(dcast(GP,ensg~app, value.var = c("Zcn")),dcast(GP,ensg~app, value.var = c("Zjs"))[,2:7],dcast(GP,ensg~app, value.var = c("Zsb"))[,2:7])
  names(GG)[2:ncol(GG)]=c(paste0(c("ACAT","GBJ","GHC","minP","TSQ","VCT"),"mc"),paste0(c("ACAT","GBJ","GHC","minP","TSQ","VCT"),"mj"),
                          paste0(c("ACAT","GBJ","GHC","minP","TSQ","VCT"),"ms"))
  return(GG)
}

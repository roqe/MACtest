#' @export
decorAPP=function(HA,PS,mc){
  # decorrelation approach
  MAp=lapply(PS, function(r){ if(all(is.na(r$P_hat))) return(NA); return(r$P_hat$a) })
  MBp=lapply(PS, function(r){ if(all(is.na(r$P_hat))) return(NA); return(r$P_hat$b) })
  CPp=ppp(data.table(za=unlist(MAp),zb=unlist(MBp)),mode = "z")
  ii=Reduce(sum,c(1,sapply(MAp,length)),accumulate = T)
  CSp=parallel::mclapply(1:length(HA$Y),function(c){
    if(is.na(PS[[c]]$P_hat)) return(NULL); return(LCN(c,CPp,PS[[c]]$P_hat$R,ii))
  },mc.cores = mc,mc.preschedule = T,mc.cleanup = T)
  # CSp=foreach::foreach(c=1:length(HA$Y)) %dopar% { LCN(c,CPp,PS[[c]]$P_hat$R,ii) }
  names(CSp)=names(HA$Y)
  LQp=data.table::rbindlist(CSp,idcol = "ensg")
  names(LQp)=c("ensg",paste0("BJ",c("pc","pj","ps")),paste0("HC",c("pc","pj","ps")),paste0("MP",c("pc","pj","ps")),paste0("AC",c("pc","pj","ps")))

  # decorrelation approach, skip decorrelation
  # MAm=sapply(PS, function(r){ if(all(is.na(r$M_hat))) return(NA); return(r$M_hat$a) })
  # MBm=sapply(PS, function(r){ if(all(is.na(r$M_hat))) return(NA); return(r$M_hat$b) })
  # CPm=ppp(list(za=unlist(MAm),zb=unlist(MBm)),mode = "z")
  # ii=Reduce(sum,c(1,sapply(MAm,length)),accumulate = T)
  # CSm=parallel::mclapply(1:length(HA$Y),function(c){
  #   if(is.na(PS[[c]]$M_hat)) return(NULL); return(LCN(c,CPm,PS[[c]]$M_hat$R,ii))
  # },mc.cores = mc,mc.preschedule = T,mc.cleanup = T)
  # # CSm=foreach::foreach(c=1:length(HA$Y)) %dopar% { LCN(c,CPm,PS[[c]]$M_hat$R,ii) }
  # names(CSm)=names(HA$Y)
  # LQm=data.table::rbindlist(CSm,idcol = "ensg")
  # names(LQm)=c("ensg",paste0("BJ",c("mc","ms","mj")),paste0("HC",c("mc","ms","mj")),paste0("MP",c("mc","ms","mj")),paste0("AC",c("mc","ms","mj")))
  # # print(Sys.time()-t1)

  return(LQp)
}

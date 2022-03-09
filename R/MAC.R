#' @param HA dataset
#' @param mc Number of cores for parallel computing, default=5.
#' @param prec The unit to construct empirical normal product pdf using in composite test, default=1e-3.
#' @importFrom parallel mclapply
#' @import GBJ
#' @import data.table
#' @export
#' @examples
#' HA=sim_mediation_data(hypo="HA",mm=0.1,vv=0.1,sm=10)
#' annoR=MAC(H00,mc=mc,prec)

MAC=function(HA,prec=1e-3,mc=5,fdr=T,sig_cut=0.05,la="P"){
  # pre-fitting
  PS=parallel::mclapply(1:length(HA$Y),function(c){
    return(pre_stat(HA$S,HA$M[[c]],HA$Y[[c]],as.data.frame(HA$X)))
  },mc.cores = mc,mc.preschedule = T,mc.cleanup = T)
  # PS=foreach::foreach(c=1:length(HA$Y)) %dopar% { pre_stat(HA$S,HA$M[[c]],HA$Y[[c]],HA$X) }

  # global approach
  GS=parallel::mclapply(1:length(HA$Y),function(c){
    return(GCN(HA$S,HA$M[[c]],HA$Y[[c]],as.data.frame(HA$X),PS[[c]]))
  },mc.cores = mc,mc.preschedule = T,mc.cleanup = T)
  # GS=foreach::foreach(c=1:length(HA$Y)) %dopar% { GCN(HA$S,HA$M[[c]],HA$Y[[c]],HA$X,PS[[c]]) }
  Q=lapply(GS,function(gs){ return(data.table::rbindlist(gs,idcol = "app")) })
  GQ=data.table::rbindlist(Q,idcol = "ensg")
  GP=lapply(split(GQ,GQ$app),ppp,mode = "p",prec,mc)

  # local approach, decorrelation
  MAp=lapply(PS, function(r){ if(is.na(r$P_hat)) return(NA); return(r$P_hat$a) })
  MBp=lapply(PS, function(r){ if(is.na(r$P_hat)) return(NA); return(r$P_hat$b) })
  CPp=ppp(list(za=unlist(MAp),zb=unlist(MBp)),mode = "z",prec,mc=mc)
  ii=Reduce(sum,c(1,sapply(MAp,length)),accumulate = T)
  CSp=parallel::mclapply(1:length(HA$Y),function(c){
    if(is.na(PS[[c]]$P_hat)) return(NULL); return(LCN(c,CPp,PS[[c]]$P_hat$R,ii))
  },mc.cores = mc,mc.preschedule = T,mc.cleanup = T)
  # CSp=foreach::foreach(c=1:length(HA$Y)) %dopar% { LCN(c,CPp,PS[[c]]$P_hat$R,ii) }
  LQp=data.table::rbindlist(CSp,idcol = "ensg")

  # local approach, skip decorrelation
  MAm=lapply(PS, function(r){ if(is.na(r$M_hat)) return(NA); return(r$M_hat$a) })
  MBm=lapply(PS, function(r){ if(is.na(r$M_hat)) return(NA); return(r$M_hat$b) })
  CPm=ppp(list(za=unlist(MAm),zb=unlist(MBm)),mode = "z",prec,mc=mc)
  ii=Reduce(sum,c(1,sapply(MAm,length)),accumulate = T)
  CSm=parallel::mclapply(1:length(HA$Y),function(c){
    if(is.na(PS[[c]]$M_hat)) return(NULL); return(LCN(c,CPm,PS[[c]]$M_hat$R,ii))
  },mc.cores = mc,mc.preschedule = T,mc.cleanup = T)
  # CSm=foreach::foreach(c=1:length(HA$Y)) %dopar% { LCN(c,CPm,PS[[c]]$M_hat$R,ii) }
  LQm=data.table::rbindlist(CSm,idcol = "ensg")

  # combine and annotate
  QQ=cbind(LQm,LQp[,-1],GP$VCT[,c("Pcn","Psb","Pjs")],GP$TSQ[,c("Pcn","Psb","Pjs")],
           GP$GBJ[,c("Pcn","Psb","Pjs")],GP$GHC[,c("Pcn","Psb","Pjs")],GP$mnP[,c("Pcn","Psb","Pjs")])

  names(QQ)[2:34]=c(paste0("BJ",c("mc","ms","mj")),paste0("HC",c("mc","ms","mj")),paste0("MP",c("mc","ms","mj")),
                    paste0("BJ",c("pc","ps","pj")),paste0("HC",c("pc","ps","pj")),paste0("MP",c("pc","ps","pj")),
                    paste0("VCT",c("mc","ms","mj")),paste0("TSQ",c("mc","ms","mj")),
                    paste0("GBJ",c("mc","ms","mj")),paste0("GHC",c("mc","ms","mj")),paste0("mnP",c("mc","ms","mj")))
  if(fdr){ QQ=cbind(ensg=QQ$ensg,data.table(apply(QQ[,2:34],2,p.adjust,method="fdr"))) }
  if(la=="P"){
    QQ$localSig=ifelse(QQ$BJpc<sig_cut|QQ$HCpc<sig_cut|QQ$MPpc<sig_cut,1,0)
  }else{
    QQ$localSig=ifelse(QQ$BJmc<sig_cut|QQ$HCmc<sig_cut|QQ$MPmc<sig_cut,1,0)
  }
  QQ$globalSig=ifelse(QQ$VCTmc<sig_cut|QQ$TSQmc<sig_cut|QQ$GBJmc<sig_cut|QQ$GHCmc<sig_cut|QQ$mnPmc<sig_cut,2,0)
  QQ$anno=factor(QQ$localSig+QQ$globalSig,levels = 0:3,labels = c("insignificant","consistent","diverse","mixed"))
  QQ$ensg=names(HA$Y)

  #return(QQ)
  return(list(PS=PS,GP=GP,LQp=LQp,LQm=LQm,PV=QQ))
}

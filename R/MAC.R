#' Run the main analysis.
#'
#' @param HA The dataset, can be obtained by running sim_mediation_data.
#' @param multi Apply the multivariate approach or not, default=true.
#' @param decor Apply the decorrelation approach or not, default=true.
#' @param hybrd Apply the hybrid approach or not, default=true.
#' @param mc Number of cores for parallel computing, default=5.
#' @param fdr Apply FDR correction in the final result or not, default=false.
#' @param sig_cut The significance level (alpha), default=0.1.
#' @importFrom parallel mclapply
#' @import GBJ
#' @import ACAT
#' @import data.table
#' @export
#' @examples
#' HA=sim_mediation_data(hypo="HA",mm=0.1,vv=0.1,sm=10)
#' annoR=MAC(HA)

MAC=function(HA,multi=T,decor=T,hybrd=T,mc=5,fdr=F,sig_cut=0.1){
  # pre-fitting
  # t1=Sys.time()
  PS=parallel::mclapply(1:length(HA$Y),function(c){
    return(pre_stat(HA$S,as.matrix(HA$M[[c]]),HA$Y[[c]],as.data.frame(HA$X),eiR = 1))
  },mc.cores = mc,mc.preschedule = T,mc.cleanup = T)
  # PS=foreach::foreach(c=1:length(HA$Y)) %dopar% { pre_stat(HA$S,HA$M[[c]],HA$Y[[c]],HA$X) }
  # print(Sys.time()-t1)

  if(multi){ MR=multiAPP(HA,PS,mc) }
  if(decor){ DR=decorAPP(HA,PS,mc) }
  QQ=merge(MR,DR,by = "ensg",all = T)

  if(fdr){ QQ=cbind(ensg=QQ$ensg,data.table(apply(QQ[,2:ncol(QQ)],2,p.adjust,method="fdr"))) }
  QQ$decorSig=ifelse(QQ$BJpc<sig_cut|QQ$HCpc<sig_cut|QQ$MPpc<sig_cut|QQ$ACpc<sig_cut,1,0)
  QQ$multiSig=ifelse(QQ$VCTmc<sig_cut|QQ$TSQmc<sig_cut|QQ$GBJmc<sig_cut|QQ$GHCmc<sig_cut|QQ$minPmc<sig_cut|QQ$ACATmc<sig_cut,2,0)
  QQ$anno=factor(QQ$decorSig+QQ$multiSig,levels = 0:3,labels = c("insignificant","consistent","diverse","mixed"))

  return(list(MR=MR,DR=DR,PV=QQ))
}

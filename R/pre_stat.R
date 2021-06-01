#' Pre-fitting the data, produce alpha hat and beta hat, required for every tests excepts VCT.
#'
#' @param S,M,Y,X Input data: exposure(vector), mediators(list of vectors), outcomes(list), and covariates(vectors).
#' @param mc Number of cores for parallel computing, default=5.
#' @importFrom parallel mclapply
#' @export
#' @examples
#' HA=sim_mediation_data(hypo="HA",mm=0.1,vv=0.1,sm=10)
#' PS=pre_stat(HA$S,HA$M,HA$Y,HA$X)

pre_stat=function(S,M,Y,X,mc=5){
  Q=mclapply(1:length(M),function(c){
    if(length(M)==length(Y)){
      Yc=Y[[c]]
    }else if(length(S)==length(Y)){
      Yc=Y
    }else{
      return(print("Data format error, Y is either a list of equal length of M, or a vector of equal length of S.") )
    }
    et=summary(lm(Yc~X+S+M[[c]]))
    if(sum(is.na(et$coefficients[,"Std. Error"]))!=0){
      return(list(M_hat=NA,P_hat=NA,signAB=NA))
    }else{
      M_hat=pre_fit(S, M[[c]], Yc, X)
      P=M[[c]]%*%M_hat$MPsvd$u
      P_hat=pre_fit(S, P, Yc, X)
      return(list(M_hat=M_hat,P_hat=P_hat,signAB=list(as=M_hat$as,bs=M_hat$bs)))
    }
  },mc.cores = mc,mc.cleanup = T,mc.preschedule = T)
  if(length(Y)==length(M)){ names(Q)=names(Y) }
  return(Q)
}

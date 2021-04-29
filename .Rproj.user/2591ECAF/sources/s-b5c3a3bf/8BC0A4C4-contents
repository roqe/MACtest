#' Intergrate alpha hat and beta hat into Z-statistics to apply the local approaches later
#'
#' @param PS The pre-fitting results (optional when using VCT, use for determining the signs of the statistics if available)
#' @param trans An indicator to apply decorrelation or not, default=TRUE.
#' @param mc Number of cores for parallel computing, default=5.
#' @export
#' @examples
#' HA=sim_mediation_data(hypo="HA",mm=0.1,vv=0.1,sm=10)
#' PS=pre_stat(HA$S,HA$M,HA$Y,HA$X)
#' CSp=cn_stat(PS,trans = T,mc)
#' CSm=cn_stat(PS,trans = F,mc)

cn_stat=function(PS,trans=T,mc=5){
  MA=sapply(PS, function(r){ return(ifelse(trans,r$P_hat["a"],r$M_hat["a"])) })
  MB=sapply(PS, function(r){ return(ifelse(trans,r$P_hat["b"],r$M_hat["b"])) })
  MR=sapply(PS, function(r){ return(ifelse(trans,r$P_hat["R"],r$M_hat["R"])) })
  ii=Reduce(sum,c(1,sapply(MA,length)),accumulate = T)
  Ma=unlist(MA)
  Mb=unlist(MB)
  Mcn=CompNull(Ma, Mb, mc=mc)
  Msb=Sobel(Ma, Mb)
  Mjs=JointSig(Ma, Mb)
  CS=lapply(1:length(PS),function(c){
    return(list(MR=MR[[c]],
                Mcn=Mcn$zz[ii[c]:(ii[c+1]-1)],
                Msb=Msb$zz[ii[c]:(ii[c+1]-1)],
                Mjs=Mjs$zz[ii[c]:(ii[c+1]-1)]))
  })
  names(CS)=names(PS)
  return(CS)
}

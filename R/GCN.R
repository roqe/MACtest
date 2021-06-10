#' Global approach: TSQ, GBJ, GHC, minP
#'
#' @param PS The pre-fitting results
#' @param method Other approcahes beside VCT, here we provide "TSQ","GBJ","GHC", and "minP", default="minP".
#' @param mc Number of cores for parallel computing, default=5.
#' @param prec The unit to construct empirical normal product pdf using in composite test, default=1e-3.
#' @importFrom parallel mclapply
#' @import GBJ
#' @import data.table
#' @export
#' @examples
#' HA=sim_mediation_data(hypo="HA",mm=0.1,vv=0.1,sm=10)
#' PS=pre_stat(HA$S,HA$M,HA$Y,HA$X)
#' TSQg=GCN(PS,"TSQ")
#' GBJg=GCN(PS,"GBJ")
#' GHCg=GCN(PS,"GHC")
#' mnPg=GCN(PS,"minP")

## method=c("TSQ","GBJ","GHC","minP")
GCN=function(PS,method="minP",mc=5,prec=1e-3){
  GS=parallel::mclapply(PS,function(ps){
    if(is.na(ps$M_hat)){ return(list(pa=NA,pb=NA,za=NA,zb=NA)) }
    if(method=="TSQ"){
      pa=TSQ(ps$M_hat$a,round(cov2cor(ps$M_hat$covA),digits = 2))$pv
      pb=TSQ(ps$M_hat$b,round(cov2cor(ps$M_hat$covB),digits = 2))$pv
    }else if(method=="GBJ"){
      pa=GBJ::GBJ(ps$M_hat$a,round(cov2cor(ps$M_hat$covA),digits = 2))$GBJ_pvalue
      pb=GBJ::GBJ(ps$M_hat$b,round(cov2cor(ps$M_hat$covB),digits = 2))$GBJ_pvalue
    }else if(method=="GHC"){
      pa=GBJ::GHC(ps$M_hat$a,round(cov2cor(ps$M_hat$covA),digits = 2))$GHC_pvalue
      pb=GBJ::GHC(ps$M_hat$b,round(cov2cor(ps$M_hat$covB),digits = 2))$GHC_pvalue
    }else{
      pa=GBJ::minP(ps$M_hat$a,round(cov2cor(ps$M_hat$covA),digits = 2))$minP_pvalue
      pb=GBJ::minP(ps$M_hat$b,round(cov2cor(ps$M_hat$covB),digits = 2))$minP_pvalue
    }
    return(list(pa=pa,pb=pb,za=safe_z(pa)*ps$signAB$as,zb=safe_z(pb)*ps$signAB$bs))
  },mc.cores = mc,mc.preschedule = T,mc.cleanup = T)
  return(ppp(GS,mc,prec))
}

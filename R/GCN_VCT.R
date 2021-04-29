#' Global approach: VCT
#'
#' @param S,M,Y,X Input data: exposure(vector), mediators(list of vectors), outcomes(list), and covariates(vectors).
#' @param PS The pre-fitting results (optional when using VCT, use for determining the signs of the statistics if available)
#' @param mc Number of cores for parallel computing, default=5.
#' @param prec The unit to construct empirical normal product pdf using in composite test, default=1e-3.
#' @importFrom parallel mclapply
#' @import SKAT
#' @import data.table
#' @export
#' @examples
#' HA=sim_mediation_data(hypo="HA",mm=0.1,vv=0.1,sm=10)
#' PS=pre_stat(HA$S,HA$M,HA$Y,HA$X)
#' VCTg=GCN_VCT(HA$S,HA$M,HA$Y,HA$X,PS)

GCN_VCT<-function(S, M, Y, X, PS=NULL, mc=5, prec=1e-3){
  GS=mclapply(1:length(Y),function(c){
    M.res<-resid(lm(M[[c]]~X))
    if(all(S %in% c(0,1))){
      S.res<-resid(glm(S~X,family = "binomial"))
    }else{
      S.res<-resid(lm(S~X))
    }
    pa<-TEGS(Y1=t(M.res), X=S.res, type=3, n.perm=100, factor.adaptive=TRUE, index=1:ncol(M[[c]]))$p[2]
    if (all(Y[[c]] %in% 0:1)) obj<-SKAT_Null_Model(Y[[c]]~-1+X+S, out_type="D")
    if (!all(Y[[c]] %in% 0:1)) obj<-SKAT_Null_Model(Y[[c]]~-1+X+S, out_type="C")
    pb<-SKAT(M[[c]], obj, is_check_genotype=FALSE, kernel="linear")$p.value
    if(is.null(PS)|is.na(PS[[c]]$signAB)){
      return(list(pa=pa,pb=pb,za=safe_z(pa),zb=safe_z(pb)))
    }else{
      return(list(pa=pa,pb=pb,za=safe_z(pa)*PS[[c]]$signAB$as,zb=safe_z(pb)*PS[[c]]$signAB$bs))
    }
  },mc.cores = mc,mc.preschedule = T,mc.cleanup = T)
  names(GS)=names(Y)
  return(ppp(GS,mc,prec))
}

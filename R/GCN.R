#' @export
GCN=function(S,M,Y,X,PS){
  gVCT=VCT(S,M,Y,X)
  if(is.na(PS$M_hat)){
    fit.m=lm(M~S+.,data = X)
    as=ifelse(sum(fit.m$coefficients["S",])>0,1,-1)
    fit.y=glmnet::glmnet(cbind(M,S,X), Y)
    bs=ifelse(sum(fit.y$beta[1:ncol(M),ncol(fit.y$beta)])>0,1,-1)
    return(list(VCT=list(za=safe_z(gVCT$pa)*as,zb=safe_z(gVCT$pb)*bs),
                TSQ=list(za=NA,zb=NA),
                GBJ=list(za=NA,zb=NA),
                GHC=list(za=NA,zb=NA),
                mnP=list(za=NA,zb=NA)))
  }
  MA=round(cov2cor(PS$M_hat$covA),digits = 4)
  MB=round(cov2cor(PS$M_hat$covB),digits = 4)
  aTSQ=TSQ(PS$M_hat$a,MA)$pv
  bTSQ=TSQ(PS$M_hat$b,MB)$pv
  aGBJ=GBJ::GBJ(PS$M_hat$a,MA)$GBJ_pvalue
  bGBJ=GBJ::GBJ(PS$M_hat$b,MB)$GBJ_pvalue
  aGHC=GBJ::GHC(PS$M_hat$a,MA)$GHC_pvalue
  bGHC=GBJ::GHC(PS$M_hat$b,MB)$GHC_pvalue
  amnP=GBJ::minP(PS$M_hat$a,MA)$minP_pvalue
  bmnP=GBJ::minP(PS$M_hat$b,MB)$minP_pvalue
  return(list(VCT=list(za=safe_z(gVCT$pa)*PS$signAB$as,zb=safe_z(gVCT$pb)*PS$signAB$bs),
              TSQ=list(za=safe_z(aTSQ)*PS$signAB$as,zb=safe_z(bTSQ)*PS$signAB$bs),
              GBJ=list(za=safe_z(aGBJ)*PS$signAB$as,zb=safe_z(bGBJ)*PS$signAB$bs),
              GHC=list(za=safe_z(aGHC)*PS$signAB$as,zb=safe_z(bGHC)*PS$signAB$bs),
              mnP=list(za=safe_z(amnP)*PS$signAB$as,zb=safe_z(bmnP)*PS$signAB$bs)))
}

#' @export
VCT<-function(S, M, Y, X, mc=5, prec=1e-3){
  M.res<-resid(lm(M~.,data = X))
  if(all(S %in% c(0,1))){
    S.res<-resid(glm(S~.,data = X,family = "binomial"))
  }else{
    S.res<-resid(lm(S~.,data = X))
  }
  pa=TEGSS(Y1=t(M.res), X=S.res, type=2, factor.adaptive=TRUE)$p_davs
  # if (all(S %in% 0:1)) obj<-SKAT::SKAT_Null_Model(S~-1+.,data = X, out_type="D")
  # if (!all(S %in% 0:1)) obj<-SKAT::SKAT_Null_Model(S~-1+.,data = X, out_type="C")
  # pa<-SKAT::SKAT(M, obj, is_check_genotype=FALSE, kernel="linear")$p.value
  if (all(Y %in% 0:1)) obj<-SKAT::SKAT_Null_Model(Y~-1+.+S,data = X, out_type="D")
  if (!all(Y %in% 0:1)) obj<-SKAT::SKAT_Null_Model(Y~-1+.+S,data = X, out_type="C")
  pb<-SKAT::SKAT(M, obj, is_check_genotype=FALSE, kernel="linear")$p.value
  return(list(pa=pa,pb=pb))
}

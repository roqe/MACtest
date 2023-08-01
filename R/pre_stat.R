#' @export
pre_stat=function(S,M,Y,X,eiR=0.7){
  et=lm(Y~.+S+M,data = X)
  if(sum(is.na(et$coefficients))!=0|
     sum(is.na(summary(et)$coefficients))!=0){
    return(list(M_hat=NA,P_hat=NA,signAB=NA))
  }else{
    M_hat=pre_fit(S, M, Y, X)
    ci=which(cumsum(M_hat$MPsvd$d)/sum(M_hat$MPsvd$d)<=eiR)
    while(length(ci)<2){
      eiR=eiR+0.05
      ci=which(cumsum(M_hat$MPsvd$d)/sum(M_hat$MPsvd$d)<=eiR)
    }
    P<-M%*%M_hat$MPsvd$u[,ci]
    P_hat=pre_fit(S, P, Y, X)
    return(list(M_hat=M_hat,P_hat=P_hat,signAB=list(as=M_hat$as,bs=M_hat$bs)))
  }
}

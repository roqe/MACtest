#' @export
TSQ=function(Z,cm){
  TS=as.numeric(t(Z)%*%MASS::ginv(cm,tol=1e-8)%*%Z)
  r=Matrix::rankMatrix(cm,tol=1e-8)[1]
  pv=pchisq(TS,r,lower.tail = F)
  return(list(TSQ=TS,pv=pv))
}

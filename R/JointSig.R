#' @export
JointSig=function(a, b){
  zz=apply(cbind(abs(a),abs(b)),1,min)*ifelse(a*b>0,1,-1)
  pp=2*pnorm(-abs(zz))
  return(list(pp=pp,zz=zz))
}

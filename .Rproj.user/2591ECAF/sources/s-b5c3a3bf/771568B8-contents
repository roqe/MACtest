TSQ=function(Z,cm){
  TS=t(Z)%*%solve(cm)%*%Z
  pv=pchisq(TS,nrow(cm),lower.tail = F)
  return(list(TS=TS,pv=pv))
}

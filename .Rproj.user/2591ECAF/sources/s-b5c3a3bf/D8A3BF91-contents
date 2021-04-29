para=function(hypo,sm,mm,vv,m){
  if(hypo=="H00"){
    alpha=beta=c(rep(0,m))
  }else if(hypo=="G0A"){
    alpha=c(rnorm(sm,mm,vv),rep(0,(m-sm)))
    beta=c(rep(0,m))
  }else if(hypo=="G0B"){
    alpha=c(rep(0,m))
    beta=c(rnorm(sm,mm,vv),rep(0,(m-sm)))
  }else if(hypo=="HA"){
    alpha=c(rnorm(sm,mm,vv),rep(0,(m-sm)))# ##equals to the length of m
    beta=c(rnorm(sm,mm,vv),rep(0,(m-sm)))
  }else{
    alpha=c(rep(0,(m-sm)),rnorm(sm,mm,vv))# ##equals to the length of m
    beta=c(rnorm(sm,mm,vv),rep(0,(m-sm)))
  }
  return(list(a=alpha,b=beta))
}

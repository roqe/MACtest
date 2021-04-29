rr=function(pa,pb,V){
  path=pa[1]/pb[1]
  der=matrix((pa[-1]/pa[1])-(pb[-1]/pb[1]),nrow=1)
  std.error=sqrt(der%*%V%*%t(der))
  pv=2*min(pnorm(abs(log(path)/std.error)),1-pnorm(abs(log(path)/std.error)))
  return(c(path, exp(log(path)-1.96*std.error), exp(log(path)+1.96*std.error), pv))
}

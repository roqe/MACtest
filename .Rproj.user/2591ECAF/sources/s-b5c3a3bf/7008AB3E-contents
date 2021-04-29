rd=function(pa,pb,V){
  path=pa[1]-pb[1]
  der=matrix(pa[-1]-pb[-1], nrow=1)
  std.error=sqrt(der%*%V%*%t(der))
  pv=2*min(pnorm(abs(path/std.error)),1-pnorm(abs(path/std.error)))
  return(c(path, path-1.96*std.error, path+1.96*std.error, pv))
}

corrPQ=function(p,q,a.hat,a.sd,b.hat,b.sd,covA,covB){
  N=covA[p,q]*covB[p,q]+a.hat[p]*a.hat[q]*covB[p,q]+b.hat[p]*b.hat[q]*covA[p,q]
  D1=sqrt(a.sd[p]^2*b.sd[p]^2+a.hat[p]^2*b.sd[p]^2+b.hat[p]^2*a.sd[p]^2)
  D2=sqrt(a.sd[q]^2*b.sd[q]^2+a.hat[q]^2*b.sd[q]^2+b.hat[q]^2*a.sd[q]^2)
  return(N/D1/D2)
}

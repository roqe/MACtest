#' @export
CompTest<-function(a, b, prec=1e-3,mc=5){
  ab=abs(a*b)
  upb=ceiling(max(ab))*2
  B=min(upb/prec,10000)
  sq=upb*(1:B)/B
  pdf=besselK(x=sq, nu=0)
  sumpdf=sum(pdf)
  ppc=unlist(parallel::mclapply(ab/sqrt(var(a)+var(b)-1),function(cut){
    return(sum(pdf[sq>cut]))
  },mc.cores = mc,mc.cleanup = T))/sumpdf
  min_pp=min(ppc[ppc>0])
  ppc[ppc<=0]=min_pp
  ppc[ppc>=1]=0.999999
  return(list(pp=ppc,zz=safe_z(ppc)*sign(a*b)))
}

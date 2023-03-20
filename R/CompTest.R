#' @export
CompTest<-function(a, b, prec=1e-3,mc=5,qcut=0.2){
  ab=abs(a*b)
  upb=ceiling(max(ab))*2
  B=min(upb/prec,10000)
  pdf=produce_pdf(upb=upb,B=B)
  sq=upb*1:B/B
  sumpdf=sum(pdf)
  #pp0=unlist(parallel::mclapply(ab,myp,sq,pdf,mc.cores = mc,mc.cleanup = T))/sumpdf
  ppc=unlist(parallel::mclapply(ab/sqrt(sd(a)^2+sd(b)^2-1),myp,sq,pdf,mc.cores = mc,mc.cleanup = T))/sumpdf
  min_pp=min(ppc[ppc>0])
  ppc[ppc<=0]=min_pp
  ppc[ppc>=1]=0.999999
  return(list(pp=ppc,zz=safe_z(ppc)*ifelse(a*b>0,1,-1)))
}

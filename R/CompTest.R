#' @export
CompTest<-function(a, b, prec=1e-3,mc=5,rn=0,qcut=0.1){
  pp=zz=rep(NA,length(a))
  nx=(is.na(a)|is.na(b))
  a=a[!nx]
  b=b[!nx]
  if(rn!=0){
    pa=pnorm(abs(a),lower.tail = F) * 2
    pb=pnorm(abs(b),lower.tail = F) * 2
    indd=qvalue(pa)$qvalues<qcut&qvalue(pb)$qvalues<qcut
    a_null=a[!indd]
    b_null=b[!indd]
  }
  sa=unlist(ifelse(rn==0,list(a),list(c(rep(a_null,rn),a))))
  sb=unlist(ifelse(rn==0,list(b),list(c(rep(b_null,rn),b))))
  ab=abs(a*b)
  upb=ceiling(max(ab))*2
  B=min(upb/prec,10000)
  pdf=produce_pdf(upb=upb,B=B)
  sq=upb*1:B/B
  sumpdf=sum(pdf)
  # pp0=unlist(parallel::mclapply(ab,myp,sq,pdf,mc.cores = mc,mc.cleanup = T))/sumpdf
  # pp1=unlist(parallel::mclapply(ab/sd(a),myp,sq,pdf,mc.cores = mc,mc.cleanup = T))/sumpdf
  # pp2=unlist(parallel::mclapply(ab/sd(b),myp,sq,pdf,mc.cores = mc,mc.cleanup = T))/sumpdf
  # pp.comp<-pp1+pp2-pp0
  pp.comp=unlist(parallel::mclapply(ab/sqrt(sd(sa)^2+sd(sb)^2-1),myp,sq,pdf,mc.cores = mc,mc.cleanup = T))/sumpdf
  min_pp=min(pp.comp[pp.comp>0])
  pp.comp[pp.comp<=0]=min_pp
  pp.comp[pp.comp>=1]=0.999999
  pp[!nx]=pp.comp
  zz[!nx]=safe_z(pp.comp)*ifelse(a*b>0,1,-1)
  return(list(pp=pp,zz=zz))
}

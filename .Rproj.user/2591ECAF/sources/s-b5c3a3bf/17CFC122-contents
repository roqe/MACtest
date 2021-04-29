CompNull<-function(a, b, mc=5, prec=1e-3){
  pp=zz=rep(NA,length(a))
  nx=(is.na(a)|is.na(b))
  a=a[!nx]
  b=b[!nx]
  ab=abs(a*b)
  upb=ceiling(max(ab))
  B=min(upb/prec,80000)
  pdf=produce_pdf(upb=upb,B=B)
  pp0<-unlist(mclapply(ab,myp,upb,B,pdf,mc.cores = mc))
  pp1<-unlist(mclapply(ab/sd(a),myp,upb,B,pdf,mc.cores = mc))
  pp2<-unlist(mclapply(ab/sd(b),myp,upb,B,pdf,mc.cores = mc))
  pp.comp<-pp1+pp2-pp0
  min_pp=min(pp.comp[pp.comp>0])
  pp.comp[pp.comp<=0]=min_pp
  pp.comp[pp.comp>1]=1
  pp[!nx]=pp.comp
  zz[!nx]=safe_z(pp.comp)*ifelse(a*b>0,1,-1)
  return(list(pp=pp,zz=zz))
}

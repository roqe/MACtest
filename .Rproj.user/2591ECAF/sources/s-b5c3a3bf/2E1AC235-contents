# todo: add lasso if n<m
zpm=function(S,M,Y,X){
  p<-ncol(M)
  n<-nrow(M)
  if(anyNA(X)){
    print("warning: missing covariate data")
    X<-matrix(rep(1, n))
    q<-1
  } else q<-dim(X)[2]+1
  fit.m0<-lm(M~.+S, data = X)
  if (all(Y %in% 0:1)) fit.y0<-glm(Y~.+S+M, family="binomial", data = X)
  if (!all(Y %in% 0:1)) fit.y0<-lm(Y~.+S+M, data = X)
  a.hat<-fit.m0$coef["S",]
  cov.a<-vcov(fit.m0)[1:p*(q+1), 1:p*(q+1)]
  ca.svd<-svd(cov.a)
  a.sum<-sum(t(a.hat)%*%ca.svd$u%*%diag(1/sqrt(ca.svd$d),ncol=length(a.hat)))
  b.hat<-fit.y0$coef[(q+2):(q+p+1)]
  cov.b<-vcov(fit.y0)[(q+2):(q+p+1), (q+2):(q+p+1)]
  cb.svd<-svd(cov.b)
  b.sum<-sum(t(b.hat)%*%cb.svd$u%*%diag(1/sqrt(cb.svd$d),ncol=length(b.hat)))
  return(list(as=ifelse(a.sum>0,1,-1),bs=ifelse(b.sum>0,1,-1)))
}

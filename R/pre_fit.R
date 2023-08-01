#' @export
pre_fit=function(S,M,Y,X=NULL){
  if(is.null(X)){
    print("warning: missing covariate data")
    X<-data.frame(rep(1, length(S)))
  }

  fit.m<-lm(M~S+.,data = X)
  if(all(Y %in% 0:1)){
    fit.y<-glm(Y~M+S+.,data = X, family="binomial")
  }else{
    fit.y<-glm(Y~M+S+.,data = X)
  }
  fit.m.sum<-summary(fit.m)
  fit.y.sum<-summary(fit.y)

  A=sapply(fit.m.sum,function(fms){
    a_hat<-fms$coefficients["S", 1]
    a_sd<-fms$coefficients["S", 2]
    return(c(a_hat,a_sd))
  })
  a.hat<-A[1,]
  a.sd<-A[2,]
  a=a.hat/a.sd
  indA=nrow(fit.m$coefficients)*c(0:(ncol(M)-1))+2
  covA=as.matrix(vcov(fit.m)[indA,indA])

  indB=2:(ncol(M)+1)
  b.hat<-fit.y.sum$coefficients[indB, 1]
  b.sd<-fit.y.sum$coefficients[indB, 2]
  b=b.hat/b.sd
  covB=as.matrix(vcov(fit.y)[indB,indB])

  MPsvd<-svd(cov(fit.m$residuals))
  corrR=diag(ncol(M))
  for(p in 1:ncol(M)){
    for(q in 1:p){
      corrR[q,p]=corrR[p,q]=corrPQ(p,q,a.hat,a.sd,b.hat,b.sd,covA,covB)
    }
  }
  ca.svd<-svd(covA)
  a.sum<-sum(t(a.hat)%*%ca.svd$u%*%diag(1/sqrt(ca.svd$d),ncol=length(a.hat)))
  cb.svd<-svd(covB)
  b.sum<-sum(t(b.hat)%*%cb.svd$u%*%diag(1/sqrt(cb.svd$d),ncol=length(b.hat)))

  return(list(a=a,b=b,R=corrR,covA=covA,covB=covB,MPsvd=MPsvd,as=sign(a.sum),bs=sign(b.sum)))
}

pre_fit=function(S,M,Y,X){
  if(anyNA(X)){
    print("warning: missing covariate data")
    X=matrix(rep(1, length(Y)))
  }

  fit.m<-lm(M~X+S)
  if(all(Y %in% 0:1)){
    fit.y<-glm(Y~X+S+M, family="binomial")
  }else{
    fit.y<-lm(Y~X+S+M)
  }

  fit.m.sum<-summary(fit.m)
  fit.y.sum<-summary(fit.y)
  lengX=0
  for(x in 1:ncol(X)){
    lengX=lengX+1
    if(class(X[,x])=="factor"){
      lengX=lengX+(length(levels(X[,x]))-2)
    }
  }
  if(ncol(M)!=1){
    A=sapply(fit.m.sum,function(fms){
      a_hat<-fms$coefficients[lengX+2, 1]
      a_sd<-fms$coefficients[lengX+2, 2]
      return(c(a_hat,a_sd))
    })
    a.hat<-A[1,]
    a.sd<-A[2,]
  }else{
    a.hat=fit.m.sum$coefficients[lengX+2, 1]
    a.sd=fit.m.sum$coefficients[lengX+2, 2]
  }
  MPsvd<-svd(cov(fit.m$residuals))
  indB=(lengX+3):(lengX+2+ncol(M))
  b.hat<-fit.y.sum$coefficients[indB, 1]
  b.sd<-fit.y.sum$coefficients[indB, 2]
  a=a.hat/a.sd
  b=b.hat/b.sd
  indA=(lengX+2)*c(1:ncol(M))
  covA=as.matrix(vcov(fit.m)[indA,indA])
  covB=as.matrix(vcov(fit.y)[indB,indB])
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

  return(list(a=a,b=b,R=corrR,covA=covA,covB=covB,MPsvd=MPsvd,as=ifelse(a.sum>0,1,-1),bs=ifelse(b.sum>0,1,-1)))
}

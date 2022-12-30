#' @export
TEGSS<-function(
  Mres=NA, 			# p-by-n matrix containing expression values.
  Sres=NA, 			# n numeric values indicating two phenotypic status (0/1).
  type=2, 		# type of working covariance (1: independence; 2: unstructured; 3: estimated cov with factor analyses; 4: compound symmetry).
  n.perm=NA, 		# no. of permutation.
  factor.adaptive=TRUE, 	# for type=3 only; TRUE if using adaptive approach, FALSE otherwise.
  centered=FALSE
){
  n<-dim(Mres)[2]
  p<-dim(Mres)[1]

  ## Calculate (conditional) sample covariance
  if (!centered) Mbar<-apply(Mres, 1, mean)
  if (centered) Mbar<-rep(0, dim(Mres)[1])
  Mres.cent<-Mres-Mbar
  covM<-cov(t(Mres.cent))

  ## Calculate working covariance
  if(n/p<2){
    covM<-covM+diag(sort(diag(covM))[round(quantile(1:p, probs=0.05))], p)
    n.perm=100
  }

  if (type==1){
    covM.ind<-matrix(0, p, p)
    diag(covM.ind)<-diag(covM)
    v.work<-covM.ind
  }
  if (type==2){
    v.work=covM
  }
  if (type==3){
    if (factor.adaptive==TRUE) {
      sv<-svd(covM)
      n.factor<-max(2,sum(cumsum(sv$d)/sum(sv$d)<0.8)+1)
      sel<-1:n.factor
      factor2<-(sv$u[,sel])%*%diag(sqrt(sv$d[sel]))
      covM.est<-factor2%*%t(factor2)
      diag(covM.est)<-diag(covM)
    } else {
      fa<-factanal(covmat=covM, factors=2)
      cor.est<-(fa$loadings[,])%*%t(fa$loadings[,])+diag(fa$uniqueness)
      covM.est<-diag(sqrt(diag(covM)))%*%cor.est%*%diag(sqrt(diag(covM)))
    }
    v.work<-covM.est
  }
  if (type==4){
    vij<-(sum(covM)-sum(diag(covM)))/(p^2-p)
    exch<-matrix(vij, p, p)
    diag(exch)<-1
    covM.csy<-diag(sqrt(diag(covM)), p)%*%exch%*%diag(sqrt(diag(covM)), p)
    v.work<-covM.csy
  }

  ## Calculate observed Q
  Mres.bar<-apply(Mres, 1, mean)
  v.work.inv<-chol2inv(chol((v.work)))
  Q.half=apply(sapply(1:n,function(i){ (Sres[i]*t(Mres[,i]-Mres.bar)%*%v.work.inv) }),1,sum)
  Q=sum(Q.half^2)
  LD=lapply(1:n,function(i){
    return((diag(Sres[i],nrow = nrow(v.work))%*%v.work.inv%*%diag(Sres[i],nrow = nrow(v.work))))
  })
  LL=Reduce("+",LD)
  LR=eigen(LL)$values
  accu=1e-2
  p_davs=CompQuadForm::davies(Q,LR,acc=accu)$Qq
  if(p_davs>=1) n.perm=100
  while(p_davs==0){
    accu=accu*1e-2
    p_davs=CompQuadForm::davies(Q,LR,acc=accu)$Qq
    if(accu<1e-20){
      n.perm=100
      break
    }
  }
  p_perm=p_rsch=NA

  ## Calculate Q under the null by permutation
  if(!is.na(n.perm)){
    Qj.perm=sapply(1:n.perm,function(pp){
      Sres.perm<-sample(Sres)
      if (type==1){
        Q.temp.perm<-t(as.numeric(as.matrix(Mres-Mres.bar))/rep(1, n))*rep(Sres.perm, each=p)
        Q.half.perm<-apply(matrix(Q.temp.perm, nrow=p), 1, sum)
      }else{
        Q.half.perm=apply(sapply(1:n,function(i){
          (Sres.perm[i]*t(Mres[,i]-Mres.bar)%*%v.work.inv)
        }),1,sum)
      }
      return(Q.half.perm^2)
    })

    Q.perm=apply(Qj.perm,2,sum)
    ## Calculate p-value (p1: permutation p-value; p2: scaled chi-sqaure p-value)
    p_perm<-mean(Q.perm>Q)
    p_rsch<-pchisq(Q*2*mean(Q.perm)/var(Q.perm), df=2*mean(Q.perm)^2/var(Q.perm), lower.tail=F)
  }

  if (type==1) {work="Indp"; work2="Working independence"}
  if (type==2) {work="Unst"; work2="Unstructured"}
  if (type==3) {work="Un.e"; work2="V estimated by 2 factors"}
  if (type==4) {work="Exch"; work2="Compound symmetry"}

  ## Output the results
  A=list(Q=Q,stru=work2,matx=v.work,p_rsch=p_rsch,p_perm=p_perm,p_davs=p_davs)
  return(A)
}

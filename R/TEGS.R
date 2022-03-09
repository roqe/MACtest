#' @export
TEGS<-function(
  Y1=NA, 			# p-by-n matrix containing expression values.
  X=NA, 			# n numeric values indicating two phenotypic status (0/1).
  type=2, 		# type of working covariance (1: independence; 2: unstructured; 3: estimated cov with factor analyses; 4: compound symmetry).
  n.perm=100, 		# no. of permutation.
  factor.adaptive=TRUE, 	# for type=3 only; TRUE if using adaptive approach, FALSE otherwise.
  centered=FALSE,
  index=NA		# index for gene(s)
){
  n<-dim(Y1)[2]
  p<-dim(Y1)[1]
  Y<-as.numeric(as.matrix(Y1))

  ## Calculate (conditional) sample covariance
  if (!centered) ybar<-apply(Y1, 1, mean)
  if (centered) ybar<-rep(0, dim(Y1)[1])
  Y1.cent<-Y1-ybar

  covY<-cor(t(Y1.cent))	# Note: use correlation instead of covariance
  corY<-cor(t(Y1.cent))	#

  dummy<-model.matrix(~as.factor(index))
  dummy[,1]<-2-apply(dummy, 1, sum)
  select<-dummy%*%t(dummy)
  corY[select==0]<-0
  covY[select==0]<-0

  ## Calculate working covariance
  if (type==1){
    covY.ind<-matrix(0, p, p)
    diag(covY.ind)<-diag(covY)
    diag(covY.ind)<-1
    v.work<-covY.ind
  }
  if (type==2){
    covY.stb<-covY+diag(sort(diag(covY))[round(quantile(1:p, probs=0.05))], p)
    v.work<-covY.stb
  }
  if (type==3){
    covY.stb<-covY+diag(sort(diag(covY))[round(quantile(1:p, probs=0.05))], p)

    if (factor.adaptive==TRUE) {
      sv<-svd(covY)
      n.factor<-sum(cumsum(sv$d)/sum(sv$d)<0.8)+1
      sel<-1:n.factor
      factor2<-(sv$u[,sel])%*%diag(sqrt(sv$d[sel]))
      covY.est<-factor2%*%t(factor2)
      diag(covY.est)<-diag(covY.stb)
    } else {
      fa<-factanal(covmat=covY.stb, factors=2)
      cor.est<-(fa$loadings[,])%*%t(fa$loadings[,])+diag(fa$uniqueness)
      covY.est<-diag(sqrt(diag(covY.stb)))%*%cor.est%*%diag(sqrt(diag(covY.stb)))
    }

    v.work<-covY.est
  }
  if (type==4){
    covY.stb<-covY+diag(sort(diag(covY))[round(quantile(1:p, probs=0.05))], p)
    vij<-(sum(corY)-sum(diag(corY)))/(p^2-p)
    exch<-matrix(vij, p, p)
    diag(exch)<-1
    covY.csy<-diag(sqrt(diag(covY)), p)%*%exch%*%diag(sqrt(diag(covY)), p)
    v.work<-covY.csy
  }

  ## Calculate observed Q
  Y1.bar<-apply(Y1, 1, mean)
  v.work.inv<-chol2inv(chol(v.work))
  Q.half<-0
  for (i in 1:n) Q.half<-Q.half+(X[i]*t(Y1[,i]-Y1.bar)%*%v.work.inv)
  Q<-sum(Q.half^2)

  ## Calculate Q under the null by permutation
  Q.perm<-NULL
  Qj.perm<-NULL
  for (pp in 1:n.perm){

    #set.seed(1357+pp)
    x.perm<-sample(X)

    # calculate (conditional) sample covariance
    if (!centered) ybar.perm<-apply(Y1, 1, mean)
    if (centered) ybar.perm<-rep(0, dim(Y1)[1])

    Y1.cent.perm<-Y1-ybar.perm

    covY.perm<-cor(t(Y1.cent.perm))
    corY.perm<-cor(t(Y1.cent.perm))

    corY.perm[select==0]<-0
    covY.perm[select==0]<-0

    # working unstructured
    if (type==2){
      covY.stb.perm<-covY.perm+diag(sort(diag(covY.perm))[round(quantile(1:p, probs=0.05))], p)
      v.work.perm<-covY.stb.perm
    }

    # working correlation estimated by factor analysis (2)
    if (type==3){
      covY.stb.perm<-covY.perm+diag(sort(diag(covY.perm))[round(quantile(1:p, probs=0.05))], p)

      if (factor.adaptive==TRUE) {
        sv.perm<-svd(covY.perm)
        n.factor<-sum(cumsum(sv.perm$d)/sum(sv.perm$d)<0.8)+1
        sel<-1:n.factor
        factor2.perm<-(sv.perm$u[,sel])%*%diag(sqrt(sv.perm$d[sel]))
        covY.est.perm<-factor2.perm%*%t(factor2.perm)
        diag(covY.est.perm)<-diag(covY.stb.perm)
      } else {
        fa.perm<-factanal(covmat=covY.stb.perm, factors=2)
        cor.est.perm<-(fa.perm$loadings[,])%*%t(fa.perm$loadings[,])+diag(fa.perm$uniqueness)
        covY.est.perm<-diag(sqrt(diag(covY.stb.perm)))%*%cor.est.perm%*%diag(sqrt(diag(covY.stb.perm)))
      }

      v.work.perm<-covY.est.perm
    }

    # working exchangeable
    if (type==4){
      covY.stb.perm<-covY.perm+diag(sort(diag(covY.perm))[round(quantile(1:p, probs=0.05))], p)
      vij.perm<-(sum(corY.perm)-sum(diag(corY.perm)))/(p^2-p)
      exch.perm<-matrix(vij.perm, p, p)
      diag(exch.perm)<-1
      covY.csy.perm<-diag(sqrt(diag(covY.perm)), p)%*%exch.perm%*%diag(sqrt(diag(covY.perm)), p)
      v.work.perm<-covY.csy.perm
    }

    # calculate permuted Q
    if (type!=1){
      Y1.bar<-apply(Y1, 1, mean)
      v.work.perm.inv<-chol2inv(chol(v.work.perm))
      Q.half.perm<-0
      for (i in 1:n) Q.half.perm<-Q.half.perm+(x.perm[i]*t(Y1[,i]-Y1.bar)%*%v.work.perm.inv)
      Q.perm<-c(Q.perm, sum(Q.half.perm^2))
      Qj.perm<-rbind(Qj.perm, as.vector(Q.half.perm^2))
    }
    if (type==1){
      #			Q.temp.perm<-t(as.numeric(as.matrix(Y1-Y1.bar))/rep(diag(covY.perm), n))*rep(x.perm, each=p)
      Q.temp.perm<-t(as.numeric(as.matrix(Y1-Y1.bar))/rep(1, n))*rep(x.perm, each=p)
      Q.half.perm<-apply(matrix(Q.temp.perm, nrow=p), 1, sum)
      Q.perm<-c(Q.perm, sum(Q.half.perm^2))
    }
  }
  Q.perm.sub<-Q.perm

  ## Calculate p-value (p1: permutation p-value; p2: scaled chi-sqaure p-value)
  p1<-mean(Q.perm>Q)
  p2<-pchisq(Q*2*mean(Q.perm.sub)/var(Q.perm.sub), df=2*mean(Q.perm.sub)^2/var(Q.perm.sub), lower.tail=F)

  if (type==1) {work="Indp"; work2="Working independence"}
  if (type==2) {work="Unst"; work2="Unstructured"}
  if (type==3) {work="Un.e"; work2="V estimated by 2 factors"}
  if (type==4) {work="Exch"; work2="Compound symmetry"}

  ## Output the results
  A<-NULL
  A$Q.perm<-Q.perm
  A$Qj.perm<-Qj.perm
  A$Qj<-Q.half^2
  A$rchisq<-(var(Q.perm.sub)/(2*mean(Q.perm.sub)))*rchisq(n.perm, df=2*mean(Q.perm.sub)^2/var(Q.perm.sub))
  A$work2<-work2
  A$work<-work
  A$p<-c(p1, p2)
  A$v.work=v.work
  A$Q=Q
  return(A)
}

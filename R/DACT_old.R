#' @import qvalue
#' @export
DACT_old=function(a,b){
  pa.ob <- pnorm(abs(a), lower.tail = F) * 2
  pb.ob <- pnorm(abs(b), lower.tail = F) * 2
  u3.ob <- pa.ob*I(pa.ob > pb.ob) + pb.ob*I(pa.ob <= pb.ob)
  p3.ob <- u3.ob^2
  qa <- qvalue::qvalue(pa.ob)
  pa.est <- qa$pi0
  qb <- qvalue::qvalue(pb.ob)
  pb.est <- qb$pi0
  w1.est <- 1 - pb.est
  w2.est <- 1 - pa.est
  w3.est <- 1 - w1.est - w2.est
  pDAT <- pa.ob * w1.est + pb.ob * w2.est + p3.ob * w3.est
  return(list(pp=as.numeric(pDAT),zz=safe_z(pDAT)))
}

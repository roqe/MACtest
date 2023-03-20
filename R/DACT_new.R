#' @import DACT
#' @export
DACT_new = function(a,b,correction="JC"){
  pi0a = min(1 - nonnullPropEst(a,0,1),1)
  pi0b = min(1 - nonnullPropEst(b,0,1),1)
  pa <- pnorm(abs(a), lower.tail = F) * 2
  pb <- pnorm(abs(b), lower.tail = F) * 2
  p.mat = cbind(pa,pb)
  p3 = (apply(p.mat,1,max))^2
  wg1 = pi0a*(1-pi0b)
  wg2 = (1-pi0a)*pi0b
  wg3 = pi0a*pi0b
  wg.sum = wg1 + wg2 + wg3
  wg.std = c(wg1,wg2,wg3)/wg.sum
  p_dact = wg.std[1]*pa + wg.std[2]*pb + wg.std[3]*p3
  if(correction == "Efron"){
    p_dact = EfronCorrect(p_dact)
  }
  if(correction == "JC"){
    p_dact = JCCorrect(p_dact)
  }
  return(list(pp=as.numeric(p_dact),zz=safe_z(p_dact)))
}

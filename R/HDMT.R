#' @import HDMT
#' @export
HDMT=function(a,b,exact=0){
  pa=pnorm(abs(a), lower.tail = F) * 2
  pb=pnorm(abs(b), lower.tail = F) * 2
  ab=cbind(pa,pb)
  pmax <- apply(ab,1,max)
  nullprop <- HDMT::null_estimation(ab)
  pnull <- HDMT::adjust_quantile(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,nullprop$alpha1,
                           nullprop$alpha2,ab,exact=0)
  jsfdr = HDMT::fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,nullprop$alpha1,
                        nullprop$alpha2,ab,exact=0)
  return(list(pmax=pmax,pnull=pnull,jsfdr=jsfdr))
}

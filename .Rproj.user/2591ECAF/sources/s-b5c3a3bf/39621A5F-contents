#' Simulate data for analysis.
#'
#' @param exp For binary exposure, enter the exposure prevalence. For continuous exposure, enter the mean and sd of a normal distribution.
#' @param confounders Confounder values.
#' @param sample_size Number of samples.
#' @param SC A vector of parameters, in order of: b0, bw, bq, bs, a0, aw, aq, d0, dw, sigma_q, sigma_s, where
#'    Q=rnorm(sample_size,d0+dwW,sigma_q),
#'    S=rnorm(sample_size,a0+awW+aqQ,sigma_s),
#'    Y=ifelse(rnorm(sample_size,b0+bwW+bqQ+bsS,1)<0,0,1)
#' @export
#' @examples
#' para=c(rep(-0.5,9),1,1)
#' dat1=sim_mediation_data(0.5,1000,para) #binary exposure
#' apply(dat1,2,mean) #the proportion of Y should not be too skew (nearly 0 or 1)
#' dat2=sim_mediation_data(c(0,1),1000,para) #continous exposure

sim_mediation_data=function(exp,sample_size,SC){
  if(length(exp)==1){
    W=c(rep(0,(1-exp)*sample_size),rep(1,exp*sample_size))
  }else{
    W=rnorm(sample_size,mean=exp[1],sd=exp[2])
  }
  Q=rnorm(sample_size,SC[8]+SC[9]*W,SC[10])
  S=rnorm(sample_size,SC[5]+SC[6]*W+SC[7]*Q,SC[11])
  Y=ifelse(rnorm(sample_size,SC[1]+SC[2]*W+SC[3]*Q+SC[4]*S,1)<0,0,1)
  return(cbind(Y,W,Q,S))
}

#' Apply mediation analysis for non-rare binary outcome with two continuous mediators
#'
#' @param dt Input data.
#' @param confounders Confounder values.
#' @param nb Number of bootstrapping. Default is 0 (no bootstrapping applied).
#' @param intv Number of intervention, 3 or 4. Default is 3.
#' @keywords Mediation analysis, Causal inference.
#' @export
#' @examples
#' para=c(rep(-0.5,9),1,1)
#' dat1=sim_mediation_data(0.5,1000,para) #binary exposure
#' apply(dat1,2,mean) #the proportion of Y should not be too skew (nearly 0 or 1)
#' res1_3=mediation_analysis(dat1)
#' tru1_3=calc_true_value(para)
#' dat2=sim_mediation_data(c(0,1),1000,para) #continous exposure
#' res2_4=mediation_analysis(dat2, intv=4, nb=500)
#' tru2_4=calc_true_value(para, intv=4)
#' summary(eIVF) #subset eIVF data for demo
#' res3_3=mediation_analysis(eIVF, confounders=c(log(36),log(26)))
#' res3_4=mediation_analysis(eIVF, confounders=c(log(36),log(26)), intv=4)

mediation_analysis=function(dt,confounders=c(),nb=0,intv=3){
  dt=as.data.frame(dt)
  colnames(dt)[1:4]=c("Y","W","Q","S")
  y.reg=glm(Y~., family = binomial(link="probit"), data=dt)
  s.reg=lm(S~.-Y, data=dt)
  q.reg=lm(Q~.-S-Y, data=dt)
  beta.hat=summary(y.reg)$coefficients[,1]
  alpha.hat=summary(s.reg)$coefficients[,1]
  delta.hat=summary(q.reg)$coefficients[,1]
  ss.hat=sqrt(mean((s.reg$residuals)^2))
  sq.hat=sqrt(mean((q.reg$residuals)^2))
  V.matrix=create_vmatrix(y.reg,s.reg,q.reg)
  total.effect=summary(glm(Y~.-Q-S, family = binomial (link="probit"), data=dt))$coefficients[,1]
  total.rd=pnorm(sum(total.effect*c(1,1,confounders)))-pnorm(sum(total.effect*c(1,0,confounders))) #when S=1 compared to S=0
  total.rr=pnorm(sum(total.effect*c(1,1,confounders)))/pnorm(sum(total.effect*c(1,0,confounders)))
  bdnp=c("lower(a)","upper(a)","pv(a)")

  if(nb>0){
    var.boot=create_var_boot(nb, dt, confounders=confounders, intv=intv)
    total=dim(var.boot)[1]
    b.RD1=c(quantile(var.boot[,1], probs=c(0.025, 0.975)))
    b.RD2=c(quantile(var.boot[,3], probs=c(0.025, 0.975)))
    b.RD3=c(quantile(var.boot[,5], probs=c(0.025, 0.975)))
    b.RR1=c(quantile(var.boot[,2], probs=c(0.025, 0.975)))
    b.RR2=c(quantile(var.boot[,4], probs=c(0.025, 0.975)))
    b.RR3=c(quantile(var.boot[,6], probs=c(0.025, 0.975)))
    #null hypothesis is 0 no effect in path
    valRD1=sum(var.boot[,1]>0)/total
    pvalRD1=min(valRD1,(1-valRD1))*2
    valRD2=sum(var.boot[,3]>0)/total
    pvalRD2=min(valRD2,(1-valRD2))*2
    valRD3=sum(var.boot[,5]>0)/total
    pvalRD3=min(valRD3,(1-valRD3))*2
    #null hypothesis is log(1) no difference in log(relative risk). This is done in the log scale since it is more normally distributed
    valRR1=sum(log(var.boot[,2])>0)/total
    pvalRR1=min(valRR1,(1-valRR1))*2
    valRR2=sum(log(var.boot[,4])>0)/total
    pvalRR2=min(valRR2,(1-valRR2))*2
    valRR3=sum(log(var.boot[,6])>0)/total
    pvalRR3=min(valRR3,(1-valRR3))*2
    bdnp=c("lower(a)","upper(a)","pv(a)","lower(b)","upper(b)","pv(b)")
  }

  theta_hat=list(beta.hat, alpha.hat, delta.hat, sq.hat, ss.hat)
  if(intv==3){
    p000=omega(theta_hat, c(0,0,0), confounders)
    p100=omega(theta_hat, c(1,0,0), confounders) #first part of difference
    p110=omega(theta_hat, c(1,1,0), confounders) #first part of difference
    p111=omega(theta_hat, c(1,1,1), confounders) #first part of difference
    omega_values=c(p000[1],p100[1],p110[1],p111[1],total.rd,total.rr)
    names(omega_values)=c("p000","p100","p110","p111","total RD","total RR")
    RD1=rd(p100,p000,V.matrix)
    RD2=rd(p110,p100,V.matrix)
    RD3=rd(p111,p110,V.matrix)
    RR1=rr(p100,p000,V.matrix)
    RR2=rr(p110,p100,V.matrix)
    RR3=rr(p111,p110,V.matrix)
    if(nb>0){
      pse_values=c(RD1,b.RD1,pvalRD1,RR1,b.RR1,pvalRR1,RD2,b.RD2,pvalRD2,RR2,b.RR2,pvalRR2,RD3,b.RD3,pvalRD3,RR3,b.RR3,pvalRR3)
    }else{
      pse_values=c(RD1,RR1,RD2,RR2,RD3,RR3)
    }
    names(pse_values)=c("RD W>Y",bdnp,"RR W>Y",bdnp,"RD W>S>Y",bdnp,"RR W>S>Y",bdnp,"RD W>QY",bdnp,"RR W>QY",bdnp)
  }else if(intv==4){
    p0000=omega(theta_hat, c(0,0,0,0), confounders)
    p1000=omega(theta_hat, c(1,0,0,0), confounders) #first part of difference
    p1100=omega(theta_hat, c(1,1,0,0), confounders) #first part of difference
    p1110=omega(theta_hat, c(1,1,1,0), confounders) #first part of difference
    p1111=omega(theta_hat, c(1,1,1,1), confounders) #first part of difference
    omega_values=c(p0000[1],p1000[1],p1100[1],p1110[1],p1111[1],total.rd,total.rr)
    names(omega_values)=c("p0000","p1000","p1100","p1110","p1111","total RD","total RR")
    RD1=rd(p1000,p0000,V.matrix)
    RD2=rd(p1100,p1000,V.matrix)
    RD3=rd(p1110,p1100,V.matrix)
    RD4=rd(p1111,p1110,V.matrix)
    RR1=rr(p1000,p0000,V.matrix)
    RR2=rr(p1100,p1000,V.matrix)
    RR3=rr(p1110,p1100,V.matrix)
    RR4=rr(p1111,p1110,V.matrix)
    if(nb>0){
      b.RD4=c(quantile(var.boot[,7], probs=c(0.025, 0.975)))
      b.RR4=c(quantile(var.boot[,8], probs=c(0.025, 0.975)))
      valRD4=sum(var.boot[,7]>0)/total
      pvalRD4=min(valRD4,(1-valRD4))*2
      valRR4=sum(log(var.boot[,8])>0)/total
      pvalRR4=min(valRR4,(1-valRR4))*2
      pse_values=c(RD1,b.RD1,pvalRD1,RR1,b.RR1,pvalRR1,RD2,b.RD2,pvalRD2,RR2,b.RR2,pvalRR2,
                RD3,b.RD3,pvalRD3,RR3,b.RR3,pvalRR3,RD4,b.RD4,pvalRD4,RR4,b.RR4,pvalRR4)
    }else{
      pse_values=c(RD1,RR1,RD2,RR2,RD3,RR3,RD4,RR4)
    }
    names(pse_values)=c("RD W>Y",bdnp,"RR W>Y",bdnp,"RD W>S>Y",bdnp,"RR W>S>Y",bdnp,"RD W>Q>Y",bdnp,"RR W>Q>Y",bdnp,"RD W>Q>S>Y",bdnp,"RR W>Q>S>Y",bdnp)
  }
  return(c(omega_values,pse_values))
}

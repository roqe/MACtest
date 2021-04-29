create_var_boot=function(nb, dt, confounders=c(), intv=3){
  var.boot=matrix(NA, nr=nb, nc=8)
  for(boot.count in 1:nb){
    dt.boot=as.data.frame(dt[sample(1:nrow(dt), replace=TRUE),])
    y.reg=glm(Y~., family = binomial (link="probit"), data=dt.boot)
    s.reg=lm(S~.-Y, data=dt.boot)
    q.reg=lm(Q~.-Y-S, data=dt.boot)
    beta.hat=summary(y.reg)$coefficients[,1]
    alpha.hat=summary(s.reg)$coefficients[,1]
    delta.hat=summary(q.reg)$coefficients[,1]
    ss.hat=sqrt(mean((s.reg$residuals)^2))
    sq.hat=sqrt(mean((q.reg$residuals)^2))
    theta_hat=list(beta.hat, alpha.hat, delta.hat, sq.hat, ss.hat)
    if(intv==3){
      p000=omega(theta_hat, c(0,0,0), confounders)
      p100=omega(theta_hat, c(1,0,0), confounders) #first part of difference
      p110=omega(theta_hat, c(1,1,0), confounders) #first part of difference
      p111=omega(theta_hat, c(1,1,1), confounders) #first part of difference
      RD1=p100[1]-p000[1]
      RD2=p110[1]-p100[1]
      RD3=p111[1]-p110[1]
      RR1=p100[1]/p000[1]
      RR2=p110[1]/p100[1]
      RR3=p111[1]/p110[1]
      var.boot[boot.count,]=c(RD1,RR1,RD2,RR2,RD3,RR3,NA,NA)
    }else if(intv==4){
      p0000=omega(theta_hat, c(0,0,0,0), confounders)
      p1000=omega(theta_hat, c(1,0,0,0), confounders) #first part of difference
      p1100=omega(theta_hat, c(1,1,0,0), confounders) #first part of difference
      p1110=omega(theta_hat, c(1,1,1,0), confounders) #first part of difference
      p1111=omega(theta_hat, c(1,1,1,1), confounders) #first part of difference
      RD1=p1000[1]-p0000[1]
      RD2=p1100[1]-p1000[1]
      RD3=p1110[1]-p1100[1]
      RD4=p1111[1]-p1110[1]
      RR1=p1000[1]/p0000[1]
      RR2=p1100[1]/p1000[1]
      RR3=p1110[1]/p1100[1]
      RR4=p1111[1]/p1110[1]
      var.boot[boot.count,]=c(RD1,RR1,RD2,RR2,RD3,RR3,RD4,RR4)
    }
  }
  return(var.boot)
}

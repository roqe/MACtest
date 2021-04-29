safe_z=function(pp){ return(ifelse(pp<8e-324,40,qnorm(pp/2,lower.tail = F))) }

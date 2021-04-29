Sobel=function(a, b){
  zz=a*b/sqrt(a^2+b^2)
  zz[a==0&b==0]=0
  pp=2*pnorm(-abs(zz))
  return(list(pp=pp,zz=zz))
}

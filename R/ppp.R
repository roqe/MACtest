#' @export
ppp=function(Q,mode="p"){
  if(sum(abs(Q$za)==40,na.rm = T)>0){
    Q$za[which(abs(Q$za)==40)]=max(abs(Q$za[abs(Q$za)!=40]),na.rm = T)*ifelse(Q$za[which(abs(Q$za)==40)]>0,1,-1)
  }
  if(sum(abs(Q$zb)==40,na.rm = T)>0){
    Q$zb[which(abs(Q$zb)==40)]=max(abs(Q$zb[abs(Q$zb)!=40]),na.rm = T)*ifelse(Q$zb[which(abs(Q$zb)==40)]>0,1,-1)
  }
  nna=(is.na(Q$za))|(is.na(Q$zb))
  Za=replace(Q$za,is.na(Q$za), 0)
  Zb=replace(Q$zb,is.na(Q$zb), 0)
  Mcn=CompTestER(Za,Zb)
  Msb=Sobel(Za,Zb)
  Mjs=JointSig(Za,Zb)
  if(mode=="z"){
    return(cbind(Q,Zcn=replace(Mcn$zz,nna,NA),Zsb=replace(Msb$zz,nna,NA),Zjs=replace(Mjs$zz,nna,NA)))
  }else{
    return(cbind(Q,Zcn=replace(Mcn$pp,nna,NA),Zsb=replace(Msb$pp,nna,NA),Zjs=replace(Mjs$pp,nna,NA)))
  }
}

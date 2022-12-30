#' @export
ppp=function(Q,mode="p",prec=1e-3,mc=5){
  if(sum(abs(Q$za)==40,na.rm = T)>0){
    Q$za[which(abs(Q$za)==40)]=max(abs(Q$za[abs(Q$za)!=40]),na.rm = T)*ifelse(Q$za[which(abs(Q$za)==40)]>0,1,-1)
  }
  if(sum(abs(Q$zb)==40,na.rm = T)>0){
    Q$zb[which(abs(Q$zb)==40)]=max(abs(Q$zb[abs(Q$zb)!=40]),na.rm = T)*ifelse(Q$zb[which(abs(Q$zb)==40)]>0,1,-1)
  }
  # Q$za[Q$za==0]=min(Q$za[Q$za!=0],na.rm = T)
  # Q$zb[Q$zb==0]=min(Q$zb[Q$zb!=0],na.rm = T)
  Mcn=CompTest(Q$za,Q$zb,prec,mc)
  Msb=Sobel(Q$za,Q$zb)
  Mjs=JointSig(Q$za,Q$zb)
  if(mode=="z"){
    return(data.table::data.table(Za=Q$za,Zb=Q$zb,Zcn=Mcn$zz,Zsb=Msb$zz,Zjs=Mjs$zz))
  }else{
    return(data.table::data.table(Za=Q$za,Zb=Q$zb,Pcn=Mcn$pp,Psb=Msb$pp,Pjs=Mjs$pp))
  }
}

#' @export
ppp=function(Q,mode="p",prec=1e-3,mc=5){
  Mcn=CompTest(Q$za,Q$zb,prec,mc)
  Msb=Sobel(Q$za,Q$zb)
  Mjs=JointSig(Q$za,Q$zb)
  if(mode=="z"){
    return(data.table::data.table(Za=Q$za,Zb=Q$zb,Zcn=Mcn$zz,Zsb=Msb$zz,Zjs=Mjs$zz))
  }else{
    return(data.table::data.table(Za=Q$za,Zb=Q$zb,Pcn=Mcn$pp,Psb=Msb$pp,Pjs=Mjs$pp))
  }
}

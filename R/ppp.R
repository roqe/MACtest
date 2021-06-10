ppp=function(GS,mc=5,ppp=1e-3){
  Q=data.table::rbindlist(GS,idcol = "ensg")
  Mcn=CompNull(Q$za,Q$zb,mc)
  Msb=Sobel(Q$za,Q$zb)
  Mjs=JointSig(Q$za,Q$zb)
  return(data.table::data.table(ensg=Q$ensg,Za=Q$za,Zb=Q$zb,Pcn=Mcn$pp,Psb=Msb$pp,Pjs=Mjs$pp))
}

#' @import GBJ
#' @import ACAT
#' @export
LCN=function(c,CP,MR,ii){
  Zcn=CP$Zcn[ii[c]:(ii[c+1]-1)]
  Zsb=CP$Zsb[ii[c]:(ii[c+1]-1)]
  Zjs=CP$Zjs[ii[c]:(ii[c+1]-1)]
  SS=round(MR,digits = 8)
  return(list(BJcn=GBJ::BJ(Zcn,SS)$BJ_pvalue,
              BJsb=GBJ::BJ(Zsb,SS)$BJ_pvalue,
              BJjs=GBJ::BJ(Zjs,SS)$BJ_pvalue,
              HCcn=GBJ::HC(Zcn,SS)$HC_pvalue,
              HCsb=GBJ::HC(Zsb,SS)$HC_pvalue,
              HCjs=GBJ::HC(Zjs,SS)$HC_pvalue,
              MPcn=GBJ::minP(Zcn,SS)$minP_pvalue,
              MPsb=GBJ::minP(Zsb,SS)$minP_pvalue,
              MPjs=GBJ::minP(Zjs,SS)$minP_pvalue,
              ACcn=ACAT::ACAT(2*pnorm(abs(Zcn),lower.tail = F)),
              ACsb=ACAT::ACAT(2*pnorm(abs(Zsb),lower.tail = F)),
              ACjs=ACAT::ACAT(2*pnorm(abs(Zjs),lower.tail = F))
              ))
}

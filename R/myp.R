myp=function(cut,upb,B,pdf=pdf){
  select=(upb*1:B/B)>cut
  pdf.sub=pdf[select]
  pval=sum(pdf.sub)/sum(pdf)
  return(pval)
}

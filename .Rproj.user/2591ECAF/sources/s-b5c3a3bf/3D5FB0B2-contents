#' Calculate true value of PSEs.
#'
#' @param SC A vector of parameters, in order of: b0, bw, bq, bs, a0, aw, aq, d0, dw, sigma_q, sigma_s.
#' @param intv Number of intervention, 3 or 4. Default is 3.
#' @export
#' @examples
#' para=c(rep(-0.5,9),1,1)
#' tru3=calc_true_value(para)
#' tru4=calc_true_value(para, intv=4)

calc_true_value=function(SC,intv=3){
  if(intv==3){
    PSEs=c(pse(SC,w1=c(1,0,0),w2=c(0,0,0)),
           pse(SC,w1=c(1,1,0),w2=c(1,0,0)),
           pse(SC,w1=c(1,1,1),w2=c(1,1,0)))
    names(PSEs)=c("RD W>Y","RR W>Y","RD W>S>Y","RR W>S>Y","RD W>QY","RR W>QY")
  }else if(intv==4){
    PSEs=c(pse(SC,w1=c(1,0,0,0),w2=c(0,0,0,0)),
           pse(SC,w1=c(1,1,0,0),w2=c(1,0,0,0)),
           pse(SC,w1=c(1,1,1,0),w2=c(1,1,0,0)),
           pse(SC,w1=c(1,1,1,1),w2=c(1,1,1,0)))
    names(PSEs)=c("RD W>Y","RR W>Y","RD W>S>Y","RR W>S>Y","RD W>Q>Y","RR W>Q>Y","RD W>Q>S>Y","RR W>Q>S>Y")
  }
  return(PSEs)
}

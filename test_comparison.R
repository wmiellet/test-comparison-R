#Make tool for manual calculation of predictive values, concordance and their respective confidence intervals
comparison <- function(x,y){
  library(dplyr)
  x <- if_else(x == '2','0','1')
  y <- if_else(y == '2','0','1')
  levs <- sort(union(x, y))
  tab <- table(factor(x,levs),factor(y,levs))
  tab2<-tab[order(tab[,2],decreasing=T),
            order(tab[2,],decreasing=T)]
  McN <- stats::mcnemar.test(tab2[])
  {total <- (tab2[1]+tab2[2]+tab2[3]+tab2[4])
    sp <- (tab2[4]/(tab2[4]+tab2[2]))
    sn <- (tab2[1]/(tab2[1]+tab2[3]))
    ppv <- (tab2[1]/(tab2[1]+tab2[2]))
    npv <- (tab2[4]/(tab2[4]+tab2[3]))
    cc <- ((tab2[1]+tab2[4])/total)}
  pyes <- (((tab2[1]+tab2[2])/total)*((tab2[1]+tab2[3])/total))
  pno <- (((tab2[4]+tab2[2])/total)*((tab2[4]+tab2[3])/total))
  pe <- (pyes+pno)
  kappa <- ((cc-pe)/(1-pe))
  SEsp <- (sqrt(sp*(1-sp)/(tab2[2]+tab2[4])))
  SEsn <- (sqrt(sn*(1-sn)/(tab2[1]+tab2[3])))
  SEppv <- (sqrt(ppv*(1-ppv)/(tab2[1]+tab2[2])))
  SEnpv <- (sqrt(npv*(1-npv)/(tab2[3]+tab2[4])))
  sp.ll <- sp-1.96*SEsp
  sp.ul <- sp+1.96*SEsp
  sn.ll <- sn-1.96*SEsn
  sn.ul <- sn+1.96*SEsn
  ppv.ll <- ppv-1.96*SEppv
  ppv.ul <- ppv+1.96*SEppv
  npv.ll <- npv-1.96*SEnpv
  npv.ul <- npv+1.96*SEnpv
  p1 <- prop.test(total*cc, total, conf.level=0.95, correct = FALSE)
  p2 <- as.numeric(p1$conf.int)
  cat(paste("-------------",
            "Positive predictive value",round(ppv*100,1),
            round(ppv.ll*100,1),round(ppv.ul*100,1),
            "Negative predictive value",round(npv*100,1),
            round(npv.ll*100,1),round(npv.ul*100,1),
            "sensitivity",round(sn*100,1),
            round(sn.ll*100,1),round(sn.ul*100,1),
            "specificity",round(sp*100,1),
            round(sp.ll*100,1),round(sp.ul*100,1),
            "Percent agreement",round(cc*100,1),
            round(p2[1:1]*100,1),round(p2[2:2]*100,1),
            "Cohen's kappa",round(kappa,2),
            "-------------",
            "",sep="\n"))
  output <- (McN)
  print(tab2)
  return(output)}

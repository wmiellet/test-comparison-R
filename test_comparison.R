#Make tool for manual calculation of predictive values, concordance and their respective confidence intervals
comparison <- function(x,y){
  library(dplyr)
  x <- if_else(x == '2','0','1')
  y <- if_else(y == '2','0','1')
  levs <- sort(union(x, y))
  tab <- table(factor(x,levs),factor(y,levs))
  tab2<-tab[order(tab[,2],decreasing=T),
            order(tab[2,],decreasing=T)]
  tab3<-tab2[order(tab2[,1],decreasing=F),
            order(tab2[1,],decreasing=F)]
  McN <- stats::mcnemar.test(tab3[])
  {total <- (tab3[1]+tab3[2]+tab3[3]+tab3[4])
    sp <- (tab3[4]/(tab3[4]+tab3[2]))
    sn <- (tab3[1]/(tab3[1]+tab3[3]))
    ppv <- (tab3[1]/(tab3[1]+tab3[2]))
    npv <- (tab3[4]/(tab3[4]+tab3[3]))
    cc <- ((tab3[1]+tab3[4])/total)}
  pyes <- (((tab3[1]+tab3[2])/total)*((tab3[1]+tab3[3])/total))
  pno <- (((tab3[4]+tab3[2])/total)*((tab3[4]+tab3[3])/total))
  pe <- (pyes+pno)
  kappa <- ((cc-pe)/(1-pe))
  SEsp <- (sqrt(sp*(1-sp)/(tab3[2]+tab3[4])))
  SEsn <- (sqrt(sn*(1-sn)/(tab3[1]+tab3[3])))
  SEppv <- (sqrt(ppv*(1-ppv)/(tab3[1]+tab3[2])))
  SEnpv <- (sqrt(npv*(1-npv)/(tab3[3]+tab3[4])))
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
  print(tab3)
  return(output)}

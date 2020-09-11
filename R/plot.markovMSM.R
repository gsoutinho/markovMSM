#' Plot for an object of class "markovMSM"
#'
#' @description Plots for the visual inspection for checking the Markov condition through the difference 
#' between the AJ and LMAJ from a starting time of the transition probabilities for each transition.
#'
#' @param x A dataframe in long format in case of the plot with the differences between AJ and LMAJ estimators
#' or a object of class "markovMSM" with the results of the AUC globalTest function.
#' @param quantileOrder Order of the quantil used in the formula of the AUC global test.
#' @param tran Transition probability from the starting state. Example: if objet 'from' of global is equal
#' 1 then tran='pstate2' means the results of the plot correspond to the difference between AJ and LMAJ
#' estimates of transition probabilites of the transition 1->2.
#' @param axis.scale Limit of the y axis of the plots.
#' @param ... For future methods.
#' 
#' @return No value is returned.
#' 
#' @examples
#' 
#' data("ebmt4")
#' db_wide <- ebmt4
#' positions=list(c(2, 3, 5, 6), c(4, 5, 6), c(4, 5, 6), c(5, 6), c(6), c())
#' namesStates =  c("Tx", "Rec", "AE", "Rec+AE", "Rel", "Death")
#' tmat <-transMatMSM(positions, namesStates)
#' timesNames = c(NA, "rec", "ae","recae", "rel", "srv")
#' status=c(NA, "rec.s", "ae.s", "recae.s","rel.s", "srv.s")
#' trans = tmat
#' db_long<- prepMSM(data=db_wide, trans, timesNames, status)
#' res2<-globalTestAUC(db_long, db_wide, from=2, to=3, quantiles=c(.05, .10, .20, .30, 0.40), 
#'                 tmat = tmat, replicas = 20, positions=positions, namesStates=namesStates,
#'                 timesNames=timesNames,status=status)
#' 
#' rm(res)

#' res<-localTestAUC(db_long, db_wide, times=30, from=3, to=5, tmat = tmat, 
#'                 replicas = 20, namesStates=namesStates, positions=positions,
#'                 timesNames=timesNames,status=status)

#' plot(res, tran='pstate2', axis.scale=c(0,.2))
#' 
#' @author Gustavo Soutinho and Luis Meira-Machado.

plot.markovMSM<- function(x, quantileOrder, tran='pstate2', axis.scale=c(-1, 1),...){
  
  
  if(inherits(x, "markovMSM") & class(x)[1] =='localTest'){
    
    
    for(c in 1:ncol(x$AJall)){
      
      #c<-4
      
      if(names(x$AJall)[c]==tran){
        
        AJ<-x$AJall[,c]
        LMAJ<-x$LMAJall[,c]
        
      }
    }
    
    res <- cbind(x$AJall$time, AJ, LMAJ)
    
    min01<-axis.scale[1]
    max01<-axis.scale[2]
    
    y<-c(min01, max01)
    
    matplot(res[,1],res[,2:3], type="s", col=1:2, xlab="Time (days)", ylab=paste(tran,"(", x$times, ",t)", 
                                                                                 sep = ""),
            lty=1, ylim=y[1:2])
    
    legend("bottomright",legend=c("AJ","LMAJ"),text.col=1:2)
    
    
  }
  
  if(inherits(x, "markovMSM") & class(x)[1] =='globalTest'){
    
    i<-quantileOrder
    
    s<-x$quantile[i]
    
    tempos.graf<-x$DIF$tComuns[x$DIF$s== s]
    
    for(a in 1:ncol(x$DIF)){
      #a<-1
      if(names(x$DIF)[a]==tran){
        ET<-x$DIF[x$DIF$s==s,a]
      }
    }
    
    #ET==x$DIF[x$DIF$s==s,6]
    
    for(b in 1:ncol(x$DIF)){
      #b<-3
      if(names(x$ET.qi2All)[b]==tran){
        
        ET.LI<-x$ET.qi2All[x$ET.qi2All$s==s,b]
        ET.LS<-x$ET.qs2All[x$ET.qs2All$s==s,b]
        
      }
    }
    
    xt <-tempos.graf
    F01 <-ET
    L01 <-ET.LI
    U01 <-ET.LS
    
    df01 <- data.frame(xt, F01, L01, U01)
    
    
    xx<-seq(from=s, to = (xt[length(xt)-2]), by = (xt[length(xt)-2]-s)/5)
    
    min<-axis.scale[1]
    max<-axis.scale[2]
    
    y<-c(min, max)
    
    plot(df01$xt, df01$F01, type = "s", ylim=y[1:2],
         xlab="Time (days)",
         ylab=paste("To=", tran,"(",s,",t): AJ-LM", sep=""),
         xaxt="n")
    axis(1,xx)
    polygon(c(df01$xt,rev(df01$x)),c(df01$L01,rev(df01$U01)),col = "grey75", border = FALSE)
    lines(df01$xt, df01$F01, lwd = 2)
    lines(df01$xt, df01$U01, col="blue",lty=2)
    lines(df01$xt, df01$L01, col="blue",lty=2)
    abline(h=0, col=2)
    
  }
  
}
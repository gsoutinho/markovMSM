#' Plot for an object of class "markovMSM"
#'
#' @description Plots for the visual inspection for checking the Markov condition
#' through the difference between the AJ and LMAJ from a starting time of the
#' transition probabilities for each transition.
#' 
#' @param x A dataframe in long format in case of the plot with the differences
#' between AJ and LMAJ estimators or a object of class "markovMSM" with the
#' results of the global.test or local.test function.
#' @param y Parameter with the value NULL.
#' @param quantileOrder Order of the quantil used in the formula of the AUC
#' global test.
#' @param tran Transition probability from the starting state. Example: if objet
#' 'from' of global is equal 1 then tran='pstate2' means the results of the plot
#' correspond to the difference between AJ and LMAJ estimates of transition
#' probabilites of the transition 1->2.
#' @param axis.scale Limit of the y axis of the plots.
#' @param difP Type of plot representing in case of x means the results of the
#' Localtest. If difP=TRUE plot depicts the discrepancies between AJ and LMAJ
#' estimators.  If difP=FALSE plot show the AJ and the LMAJ estimates. 
#' By default difP=FALSE.
#' @param ... For future methods.
#' 
#' @return No value is returned.
#' 
#' @examples
#'  
#' library(mstate)
#' 
#' data(prothr)
#' 
#' res<-local.test(db_long=prothr, times=30, from=1, to=3, replicas=5)
#' 
#' plot(x=res, tran='pstate2', axis.scale=c(-0.25,.25), difP=TRUE)
#'
#' plot(x=res, tran='pstate2', axis.scale=c(0,.25), difP=FALSE)
#'
#' res2<-global.test(db_long=prothr, db_wide = NULL, from=1, to=3, replicas=5,
#' limit=0.90, 
#'                  quantiles=c(.05, .10, .20, .30, 0.40))
#'
#' plot(x=res, quantileOrder=3, tran='pstate2', axis.scale=c(-0.05,.15))
#' 
#' @author Gustavo Soutinho and Luis Meira-Machado.


plot.markovMSM<- function(x = object, y = NULL, quantileOrder=NULL, tran='pstate2', 
                          axis.scale=c(-1, 1), difP=FALSE,...){
  
  object <- x
  
  z<-object
  
  if(inherits(z, "markovMSM") & class(z)[1] =='localTest' & difP==TRUE){
    
    s<-z$times
    
    tempos.graf<-z$DIF$tComuns[z$DIF$s== s]
    
    for(a in 1:ncol(z$DIF)){
      #a<-1
      if(names(z$DIF)[a]==tran){
        ET<-z$DIF[z$DIF$s==s,a]
      }
    }
    
    #ET==z$DIF[z$DIF$s==s,6]
    
    for(b in 1:ncol(z$DIF)){
      #b<-3
      if(names(z$ET.qi2All)[b]==tran){
        
        ET.LI<-z$ET.qi2All[z$ET.qi2All$s==s,b]
        ET.LS<-z$ET.qs2All[z$ET.qs2All$s==s,b]
        
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
         ylab=paste("p23(",s,",t): AJ-LMAJ", sep=""),
         xaxt="n")
    axis(1,xx)
    polygon(c(df01$xt,rev(df01$x)),c(df01$L01,rev(df01$U01)),col = "grey75", border = FALSE)
    lines(df01$xt, df01$F01, lwd = 2)
    lines(df01$xt, df01$U01, col="blue",lty=2)
    lines(df01$xt, df01$L01, col="blue",lty=2)
    abline(h=0, col=2)
    
  }
  
  if(inherits(z, "markovMSM") & class(z)[1] =='localTest' & difP==FALSE){
    
    
    for(c in 1:ncol(z$AJall)){
      
      #c<-4
      
      if(names(z$AJall)[c]==tran){
        
        AJ<-z$AJall[,c]
        LMAJ<-z$LMAJall[,c]
        
      }
    }
    
    res <- cbind(z$AJall$time, AJ, LMAJ)
    
    min01<-axis.scale[1]
    max01<-axis.scale[2]
    
    y<-c(min01, max01)
    
    matplot(res[,1],res[,2:3], type="s", col=1:2, xlab="Time (days)", ylab=paste(tran,"(", z$times, ",t)", 
                                                                                 sep = ""),
            lty=1, ylim=y[1:2])
    
    legend("bottomright",legend=c("AJ","LMAJ"),text.col=1:2)
    
  }
  
  if(inherits(z, "markovMSM") & class(z)[1] =='globalTest'){
    
    i<-quantileOrder
    
    s<-z$quantile[i]
    
    tempos.graf<-z$DIF$tComuns[z$DIF$s== s]
    
    for(a in 1:ncol(z$DIF)){
      #a<-1
      if(names(z$DIF)[a]==tran){
        ET<-z$DIF[z$DIF$s==s,a]
      }
    }
    
    #ET==z$DIF[z$DIF$s==s,6]
    
    for(b in 1:ncol(z$DIF)){
      #b<-3
      if(names(z$ET.qi2All)[b]==tran){
        
        ET.LI<-z$ET.qi2All[z$ET.qi2All$s==s,b]
        ET.LS<-z$ET.qs2All[z$ET.qs2All$s==s,b]
        
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
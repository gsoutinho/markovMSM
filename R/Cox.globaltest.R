#' Cox PH global test for checking the Markov condition on Multi-state Models.
#'
#' @description This function is used to obtain a global test based on Cox 
#' Proportional hazard models for each transition in order to check the Markov
#' condition.
#' 
#' @param data A data frame in long format containing the subject \code{id};
#' \code{from} corresponding to the starting state;the receiving state, 
#' \code{to}; the transition number, \code{trans}; the  starting time of the 
#' transition given by \code{Tstart}; the stopping time of the transition, 
#' \code{Tstop}, and \code{status} for the  status variable, with 1 indicating
#' an event (transition), 0 a censoring. 
#' @param transition Number of a transition to test the effect of the history
#' of the process for the Markov assumption until this particular transition.
#' By default value is NULL. In this case all the results for  each transition
#' are computed. 
#' 
#' @return An object with a list with the following outcomes:
#' \item{p.value}{p-value of Cox global tests for each transition.} 
#' \item{transition}{Number of a transition or all possible transitions to 
#' check the global test for the Markov assumption.}
#' 
#' @examples
#' 
#' data("colonMSM")
#' db_wide<-colonMSM
#' positions<-list(c(2, 3), c(3), c())
#' namesStates =  c("Alive", "Rec",  "Death")
#' tmat <-transMatMSM(positions, namesStates)
#' timesNames = c(NA, "time1","Stime")
#' status=c(NA, "event1","event")
#' trans = tmat
#' db_long<- prepMSM(data=db_wide, trans, timesNames, status)
#' res1<-Cox.globaltest(data=db_long, transition = 3)
#' res1
#' 
#' data("ebmt4")
#' db_wide <- ebmt4
#' positions=list(c(2, 3, 5, 6), c(4, 5, 6), c(4, 5, 6),
#'                c(5, 6), c(6), c())
#'
#' namesStates = c("Tx", "Rec", "AE", "Rec+AE", "Rel", "Death")
#' tmat <-transMatMSM(positions, namesStates)
#' timesNames = c(NA, "rec", "ae","recae", "rel", "srv")
#' status=c(NA, "rec.s", "ae.s", "recae.s","rel.s", "srv.s")
#' trans = tmat
#' db_long<- prepMSM(data=db_wide, trans, timesNames, status)
#' db_long$trans<-as.factor(db_long$trans)
#'
#' res2<-Cox.globaltest(data=db_long,transition = NULL)
#' res2$p.value
#' res2$transition
#'
#' res3<-Cox.globaltest(data=db_long, transition = 13)
#' res3
#'
#' @author Gustavo Soutinho and Luis Meira-Machado.


Cox.globaltest<-function(data, transition=NULL){
  
  db_long2<-data[data$Tstart!=0,]
  
  trans<-unique(db_long2$trans)
  
  if(is.null(transition)){
    
    dim<-length(trans)
    
    p.value<-rep(0,dim)
    
    
    for(i in 1:dim){
      
      c0 <-  suppressWarnings(coxph(Surv(Tstart, Tstop, status) ~ Tstart + 
             strata(trans), data = data[data$trans==trans[i],]))
      
      p.value[i]<-summary(c0)$coefficients[5]
      
    }
    
    res<-list(p.value = p.value,transition = trans)
   
    return(invisible(res))
    
  }else{
    
    c0 <- suppressWarnings(coxph(Surv(Tstart, Tstop, status) ~ Tstart+ 
          strata(trans), data = data[data$trans==transition,]))
    
    
    res <-list(p.value = summary(c0)$coefficients[5],transition = transition)

    return(invisible(res))
    
  }
  
}
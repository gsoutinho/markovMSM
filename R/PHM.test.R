#' Cox PH global test for checking the Markov condition on Multi-state 
#' Models.
#' @description This function is used to obtain a global test to check the 
#' Markov condition for each transition based on Cox Proportional hazard models.
#' @param data A data frame in long format containing the subject \code{id};
#' \code{from} corresponding to the starting state;the receiving state, 
#' \code{to}; the transition number, \code{trans}; the  starting time of the 
#' transition given by \code{Tstart}; the stopping time of the transition, 
#' \code{Tstop}, and \code{status} for the  status variable, with 1 indicating
#' an event (transition), 0 a censoring. 
#' @param from The starting state of the transition to check the Markov 
#' condition.
#' @param to The last state of the considered transition to check the 
#' Markov condition.
#' @return An object with a list with the following outcomes:
#' \item{p.value}{p-value of Cox global tests for each transition.} 
#' \item{from}{The starting state of the transition to check the Markov 
#' condition.}
#' \item{to}{The last state of the considered transition to check the 
#' Markov condition.}
#' @references Kay, R (1986). A Markov model for analyzing cancer markers and 
#' disease states in survival studies. \emph{Biometrics} 42, 457-481.
#' Soutinho G, Meira-Machado L (2021). Methods for checking the Markov condition
#' in multi-state survival data. \emph{Computational Statistics}.  
#' @examples
#' data("colonMSM")
#' db_wide<-colonMSM
#' positions<-list(c(2, 3), c(3), c())
#' namesStates =  c("Alive", "Rec",  "Death")
#' tmat <-transMatMSM(positions, namesStates)
#' timesNames = c(NA, "time1","Stime")
#' status=c(NA, "event1","event")
#' trans = tmat
#' db_long<- prepMSM(data=db_wide, trans, timesNames, status)
#' res1<-PHM.test(data=db_long, from = 2, to=3)
#' res1
#' 
#' data("ebmt4")
#' db_wide <- ebmt4
#' positions=list(c(2, 3, 5, 6), c(4, 5, 6), c(4, 5, 6),
#'                c(5, 6), c(6), c())
#' namesStates = c("Tx", "Rec", "AE", "Rec+AE", "Rel", "Death")
#' tmat <-transMatMSM(positions, namesStates)
#' timesNames = c(NA, "rec", "ae","recae", "rel", "srv")
#' status=c(NA, "rec.s", "ae.s", "recae.s","rel.s", "srv.s")
#' trans = tmat
#' db_long<- prepMSM(data=db_wide, trans, timesNames, status)
#' db_long$trans<-as.factor(db_long$trans)
#' res2<-PHM.test(data=db_long, from = 5, to=6)
#' res2$p.value
#' res2$from
#' res2$to
#' @author Gustavo Soutinho and Luis Meira-Machado.
#' 
#' @export PHM.test

PHM.test<-function(data, from, to){
  
  db_long2<-data[data$Tstart!=0,]
  
  trans<-unique(db_long2[db_long2$from==from & db_long2$to==to,'trans'])
  
  c0 <- suppressWarnings(coxph(Surv(Tstart, Tstop, status) ~ Tstart+ 
                                 strata(trans), data = db_long2[db_long2$trans==trans,]))
  
  res<-list(summary(c0)$coefficients[5][[1]],p.value=summary(c0)$coefficients[5], from = from, to=to)
  
  print(res[[1]])
  
  return(c(res[2],res[3],res[4]))
  
  #return(invisible(res))
  
}
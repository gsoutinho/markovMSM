#' Pritting fits of "markovMSM" class.
#'
#' @description Returns the formula and the methods of the Markov tests in
#' Multi-state models.
#'
#' @param x A object of "markovMSM" with the results of the AUC global or local
#' tests.
#' @param ... For future methods.
#' 
#' @return The formula and the methods of the Markov tests in Multi-state models.
#' @references Soutinho G, Meira-Machado L (2021). Methods for checking the 
#' Markov condition in multi-state survival data. \emph{Computational Statistics}. 
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
#' times=365
#' res<-AUC.test(db_long, db_wide, times=times, from=2, to=3, type='local', 
#' replicas=2, tmat = tmat)
#' print(res)
#' 
#' @author Gustavo Soutinho and Luis Meira-Machado.
#' @export print.markovMSM
#' @exportS3Method print markovMSM

print.markovMSM<- function(x, ...){
  
  #x<-res
  
  if(inherits(x, "markovMSM")){
    
    cat("Call:\n")
    print(x$call)
    
    
    cat("\nMethod:\n")
    
    if(class(x)[1] == "globalTest") method <- "AUC method global test"
    if(class(x)[1] == "localTest") method <- "AUC method local test"
    
    print(method)
    
  }
  
}
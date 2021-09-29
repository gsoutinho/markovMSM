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
#' res<-global.test(db_long, db_wide, from=2, to=3, quantiles=c(.05, .10, .20, 
#'                 .30, 0.40), tmat = tmat, replicas = 5, positions=positions, 
#'                 namesStates=namesStates,
#'                 timesNames=timesNames,status=status)
#' print(res)

#' @author Gustavo Soutinho and Luis Meira-Machado.

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
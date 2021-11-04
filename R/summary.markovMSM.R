#' Summarizing fits of "markovMSM" class.
#'
#' @description Returns the p-values of the AUC global and Local tests. Further
#' information
#' on the test are also given.
#'
#' @param object A object of "markovMSM" with the results of the AUC global or
#' local tests.
#' 
#' @param ... For future methods.

#' @return The p-values of the AUC global and Local tests. Further information
#' on the test are also given.
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
#' res$localTest
#' summary(res)
#' 
#' @author Gustavo Soutinho and Luis Meira-Machado.
#' @export summary.markovMSM
#' @exportS3Method summary markovMSM

summary.markovMSM<- function(object, ...){


  if(inherits(object, "markovMSM")){

    cat("\n")


    if(class(object)[1] == "globalTest"){

      cat('Markov Global Test', '\n')

      method <- "AUC Global Test"

      cat('Method:', method, '\n')

      cat('from state:', object$from, '\n')

      cat('Transitions:', names(object$DIF)[(2+object$from):ncol(object$DIF)], '\n')

      cat('p-values:', object$globalTest, '\n')

      cat("\n")

      cat("Estimation of p-value through Monte carlo simulation.")

      cat("\n")

      cat("Number of replicas: ", object$replicas)

      cat("\n")

      cat("Quantiles used in the Test: ", object$quantiles)

    }


    if(class(object)[1] == "localTest"){

      cat('Markov Local Test', '\n')

      method <- "AUC Local Test"

      cat('Method:', method, '\n')

      cat('Transitions from state', object$from, 'to', object$to)

      cat('p-values of the times for checking the Markov condition:', '\n')

      print(object$localTest)

      cat("Estimation of p-value through Monte carlo simulation.")
      cat("\n")
      cat("Number of replicas: ", object$replicas)
      cat("\n")
      cat("Times used in the Test: ", object$times)

    }

  }

}

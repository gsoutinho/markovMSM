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
#'
#' @examples

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
#' summary(res)
#' 
#' @author Gustavo Soutinho and Luis Meira-Machado.


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

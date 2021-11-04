#' Log-rank based test for the validity of the Markov assumption.
#' @description  Function LR.test performs the log-rank test described in
#' Titman & Putter (2020). 
#' @param db_long Multi-state data in \code{msdata} format. Should also contain
#' (dummy codings of) the relevant covariates; no factors allowed.
#' @param times Grid of time points at which to compute the statistic.
#' @param from The starting state of the transition to check the Markov 
#' condition.
#' @param to The last state of the considered transition to check the 
#' Markov condition. 
#' @param replicas Number of wild bootstrap replications to perform.
#' @param formula Right-hand side of the formula. If NULL will fit with no
#' covariates (formula="1" will also work), offset terms can also be specified.
#' @param fn A list of summary functions to be applied to the individual zbar
#' traces (or a list of lists)
#' @param fn2 A list of summary functions to be applied to the overall
#' chi-squared trace
#' @param dist Distribution of wild bootstrap random weights, either "poisson"
#' for centred Poisson (default), or "normal" for standard normal
#' @param min_time The minimum time for calculating optimal weights
#' @param other_weights Other (than optimal) weights can be specified here
#' @return LR.test returns an object of class "markovMSM", which is a list
#' with the following items: 
#' \item{localTestLR}{p-value of AUC local tests for each times and transitions.}  
#' \item{globalTestLR}{p-value of AUC global tests for each transition}
#' \item{times}{Grid of time points at which to compute the statistic.}
#' \item{replicas}{Number of wild bootstrap replications to perform.}   
#' \item{call}{Expression of the LR.test used.} 
#' @references Titman AC, Putter H (2020). General tests of the Markov property
#' in multi-state models. \emph{Biostatistics}.
#' @examples
#' set.seed(1234)
#' \dontrun{
#' data("colonMSM")
#' positions<-list(c(2, 3), c(3), c())
#' namesStates =  c("Alive", "Rec",  "Death")
#' tmat <-transMatMSM(positions, namesStates)
#' timesNames = c(NA, "time1","Stime")
#' status=c(NA, "event1","event")
#' trans = tmat
#' db_long<- prepMSM(data=colonMSM, trans, timesNames, status)
#' res<-LR.test(db_long=db_long, times=180, from = 2, to = 3, replicas = 1000)
#' res$globalTestLR
#' 
#' times<-c(73.5, 117, 223, 392, 681)
#' res2<-LR.test(db_long=prothr, times=times, from = 2, to = 3, replicas = 1000)
#' res2$localTestLR
#' res2$globalTestLR
#' 
#' res3<-LR.test(db_long=prothr, times=times, from = 2, to = 1, replicas = 1000)
#' res3$localTestLR
#' res3$globalTestLR}
#' 
#' @author Gustavo Soutinho and Luis Meira-Machado.
#' 
#' @export LR.test

LR.test<-function(db_long=db_long, times=times, from, to, replicas = 1000, formula = NULL,
                   fn = list(function(x) mean(abs(x), na.rm = TRUE)),
                   fn2 = list(function(x) mean(x, na.rm = TRUE)),
                   min_time = 0,
                   other_weights = NULL,
                   dist = c("poisson", "normal")){
  
  transition<-unique(db_long[db_long$from==from & db_long$to==to,'trans'])
  
  MT <- MarkovTest(data=db_long, id = "id", transition = transition, grid = times,  B = replicas)
  
  #MT$obs_chisq_trace 
  #MT$qualset
  #MT$orig_stat
  #MT$orig_ch_stat 
  #MT$p_stat_wb 
  #MT$p_ch_stat_wb 
  #MT$zbar
  #1-pchisq(MT$zbar,1)
  
  res<-list(localTestLR=round(1-pchisq(MT$obs_chisq_trace,1),3),
            globalTestLR=MT$p_ch_stat_wb, times=times, replicas=replicas)

  class(res) <- c("Log-rank Tests", "markovMSM")

  res$call <- match.call()
  
  return(res)

}
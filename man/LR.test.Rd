% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LR.test.R
\name{LR.test}
\alias{LR.test}
\title{Log-rank based test for the validity of the Markov assumption.}
\usage{
LR.test(
  db_long = db_long,
  times = times,
  from,
  to,
  replicas = 1000,
  formula = NULL,
  fn = list(function(x) mean(abs(x), na.rm = TRUE)),
  fn2 = list(function(x) mean(x, na.rm = TRUE)),
  min_time = 0,
  other_weights = NULL,
  dist = c("poisson", "normal")
)
}
\arguments{
\item{db_long}{Multi-state data in \code{msdata} format. Should also contain
(dummy codings of) the relevant covariates; no factors allowed.}

\item{times}{Grid of time points at which to compute the statistic.}

\item{from}{The starting state of the transition to check the Markov 
condition.}

\item{to}{The last state of the considered transition to check the 
Markov condition.}

\item{replicas}{Number of wild bootstrap replications to perform.}

\item{formula}{Right-hand side of the formula. If NULL will fit with no
covariates (formula="1" will also work), offset terms can also be specified.}

\item{fn}{A list of summary functions to be applied to the individual zbar
traces (or a list of lists)}

\item{fn2}{A list of summary functions to be applied to the overall
chi-squared trace}

\item{min_time}{The minimum time for calculating optimal weights}

\item{other_weights}{Other (than optimal) weights can be specified here}

\item{dist}{Distribution of wild bootstrap random weights, either "poisson"
for centred Poisson (default), or "normal" for standard normal}
}
\value{
LR.test returns an object of class "markovMSM", which is a list
with the following items: 
\item{localTestLR}{p-value of AUC local tests for each times and transitions.}  
\item{globalTestLR}{p-value of AUC global tests for each transition}
\item{times}{Grid of time points at which to compute the statistic.}
\item{replicas}{Number of wild bootstrap replications to perform.}   
\item{call}{Expression of the LR.test used.}
}
\description{
Function LR.test performs the log-rank test described in
Titman & Putter (2020).
}
\examples{
set.seed(1234)
\dontrun{
data("colonMSM")
positions<-list(c(2, 3), c(3), c())
namesStates =  c("Alive", "Rec",  "Death")
tmat <-transMatMSM(positions, namesStates)
timesNames = c(NA, "time1","Stime")
status=c(NA, "event1","event")
trans = tmat
db_long<- prepMSM(data=colonMSM, trans, timesNames, status)
res<-LR.test(db_long=db_long, times=180, from = 2, to = 3, replicas = 1000)
res$globalTestLR

times<-c(73.5, 117, 223, 392, 681)
res2<-LR.test(db_long=prothr, times=times, from = 2, to = 3, replicas = 1000)
res2$localTestLR
res2$globalTestLR

res3<-LR.test(db_long=prothr, times=times, from = 2, to = 1, replicas = 1000)
res3$localTestLR
res3$globalTestLR}

}
\references{
Titman AC, Putter H (2020). General tests of the Markov property
in multi-state models. \emph{Biostatistics}.
}
\author{
Gustavo Soutinho and Luis Meira-Machado.
}

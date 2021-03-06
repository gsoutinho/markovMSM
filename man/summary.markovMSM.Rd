% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.markovMSM.R
\name{summary.markovMSM}
\alias{summary.markovMSM}
\title{Summarizing fits of "markovMSM" class.}
\usage{
\method{summary}{markovMSM}(object, ...)
}
\arguments{
\item{object}{A object of "markovMSM" with the results of the AUC global or local tests.}

\item{...}{For future methods.}
}
\value{
The p-values of the AUC global and Local tests. Further information on the test are also given.
}
\description{
Returns the p-values of the AUC global and Local tests. Further information
on the test are also given.
}
\examples{
data("ebmt4")
db_wide <- ebmt4
positions=list(c(2, 3, 5, 6), c(4, 5, 6), c(4, 5, 6), c(5, 6), c(6), c())
namesStates =  c("Tx", "Rec", "AE", "Rec+AE", "Rel", "Death")
tmat <-transMatMSM(positions, namesStates)
timesNames = c(NA, "rec", "ae","recae", "rel", "srv")
status=c(NA, "rec.s", "ae.s", "recae.s","rel.s", "srv.s")
trans = tmat
db_long<- prepMSM(data=db_wide, trans, timesNames, status)
res<-global.test(db_long, db_wide, from=2, to=3, quantiles=c(.05, .10, .20, .30, 0.40), 
                tmat = tmat, replicas = 5, positions=positions, namesStates=namesStates,
                timesNames=timesNames,status=status)
summary(res)

}
\author{
Gustavo Soutinho and Luis Meira-Machado.
}

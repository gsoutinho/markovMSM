% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transMatMSM.R
\name{transMatMSM}
\alias{transMatMSM}
\title{Define transition matrix for multi-state model.}
\usage{
transMatMSM(positions, namesStates)
}
\arguments{
\item{positions}{List of possible transitions; x[[i]] consists of a vector of
state numbers reachable from state i.}

\item{namesStates}{A character vector containing the names of either the
competing risks or the states in the multi-state model specified by the
competing risks or illness-death model. names should have the same 
length as the list x (for transMat), or either K or K+1 (for trans.comprisk),
or 3 (for trans.illdeath).}
}
\value{
A transition matrix describing the states and transitions in the
multi-state model.
}
\description{
Define transition matrices for multi-state model. Specific
functions for defining  such transition matrices are pre-defined for common
multi-state models like the competing risks model and the illness-death model.
}
\examples{
data("ebmt4")
db_wide <- ebmt4
positions<-list(c(2, 3, 5, 6), c(4, 5, 6), c(4, 5, 6), c(5, 6), c(6), c())
namesStates = c("Tx", "Rec", "AE", "Rec+AE", "Rel", "Death")
transMatMSM(positions, namesStates)
}
\author{
Gustavo Soutinho and Luis Meira-Machado.
}

#' Area Under the two curves (AUC) global test for checking the Markov condition 
#' on Multi-state Models.
#'
#' @description This function is used to obtain a global test for checking the
#' Markov condition.
#'
#' @param db_long A data frame in long format containing the subject \code{id};
#' \code{from} corresponding to the starting state;the receiving state, 
#' \code{to}; the transition number, \code{trans}; the starting time of the
#' transition given by \code{Tstart}; the stopping time of the transition, 
#' \code{Tstop}, and \code{status} for the  status variable, with 1 indicating
#' an event (transition), 0 a censoring. 
#' @param db_wide Data frame in wide format in which to interpret time, status,
#' id or keep, if appropriate.
#' @param from The starting state of the transition probabilities.
#' @param to The last receiving state considered for the estimation of the
#' transition probabilities. All  the probabilities among the first and the
#' last states are also computed. 
#' @param quantiles Quantiles used in the formula of the Global test for the
#' AUC methods. By default takes the percentiles 0.05, .10, .20, .30, 0.40. 
#' @param times Times between the minimum time and the third quartile times used
#' in the formula of the Global test for the AUC methods. Default to NULL.
#' @param tmat The transition matrix for multi-state model.
#' @param replicas Number of replicas for the Monte Carlo simulation to
#' standardization of the T-statistic given by the difference of the areas of
#' AJ and LM transition probabilities estimates.
#' @param limit Percentile of the event time used as the upper bound for the
#' computation of the AUC-based test. 
#' @param positions List of possible transitions; x[[i]] consists of a vector of
#' state numbers reachable from state i.
#' @param namesStates A character vector containing the names of either the
#' competing risks or the states in the multi-state model specified by the
#' competing risks or illness-death model. names should have the same length as
#' the list x (for transMat), or either K or K+1 (for trans.comprisk), or 3
#' (for trans.illdeath). 
#' @param timesNames Either 1) a matrix or data frame of dimension n x S
#' (n being the number of individuals  and S the number of states in the
#' multi-state model), containing the times at which the states are visited or
#' last follow-up time, or 2) a character vector of length S containing the
#' column names  indicating these times. In the latter cases, some elements of
#' time may be NA, see Details
#' @param status Either 1) a matrix or data frame of dimension n x S,
#' containing, for each of the states, event indicators taking the value 1 if
#' the state is visited or 0 if it is not (censored), or 2) a character
#' vector of length S containing the column names indicating these status
#' variables. In the latter cases, some elements of status may be NA, 
#' see Details.  
#'
#' @return An object with a list with the following outcomes:
#' \item{globalTest}{p-value of AUC global tests for each transition. These
#' values are obtained through the minimum of the means of each two contiguous
#' quantiles times of the AUC global tests.} 
#' \item{localTest}{AUC local tests of the transition probability for each times
#' and transitions.}
#' \item{quantiles}{Quantiles times used for the AUC global tests.}
#' \item{times}{Times used for the AUC global tests.}
#' \item{DIF}{Differences between the AJ and the LMAJ estimates for each
#' transition probabilites from the starting state until the receiving state
#' given by only one replica where 's' represent each of the quantile times.}  
#' \item{from}{The starting state considered for the AUC global tests.}
#' \item{to}{The last receiving state considered for the the AUC Local tests.}
#' \item{ET.qiAll}{The lower limit of the diferences between the AJ and the LMAJ
#' estimates given by the Monte Carlo simulation in each transition for each "s"
#' quantile times. \code{ET.qi2All} means the same but missing values were
#' replaced by the previous diferences of estimators.} 
#' \item{ET.qsAll}{The upper limit of the diferences between the AJ and the LMAJ
#' estimates given by the Monte Carlo simulation in each transition for each "s"
#' quantile times. \code{ET.qs2All} means the same but missing values were
#' replaced by the previous diferences of estimators.} 
#' \item{replicas}{Number of replicas for the Monte Carlo simulation.}
#' \item{limit}{Percentil of the times used in the AUC global tests.}
#' 
#' @examples
#' 
#' data("colonIDM")
#' colonIDM$event1[colonIDM$time1==colonIDM$Stime & colonIDM$event1==1]<-0
#' db_wide<-colonIDM
#'
#' positions<-list(c(2, 3), c(3), c())
#' namesStates =  c("Alive", "Rec",  "Death")
#' tmat <-transMatMSM(positions, namesStates)
#' timesNames = c(NA, "time1","Stime")
#' status=c(NA, "event1","event")
#' trans = tmat
#' db_long<- prepMSM(data=db_wide, trans, timesNames, status)
#' res1<-global.test(db_long, db_wide, from=1, to=3, replicas = 5, tmat=tmat)
#' round(res1$globalTest,3)
#' round(res1$localTests,3)
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
#' 
#' res2<-global.test(db_long, db_wide, from=2, to=3, quantiles=c(.05, 0.10, .20,
#' .30, 0.40), times=NULL, tmat = tmat, replicas = 10, positions=positions, 
#'           namesStates=namesStates, timesNames=timesNames, status=status)
#'                 
#' round(res2$globalTest,3)
#' 
#' 
#' times<-c(161.1, 365.4, 688.9, 1125.8,72.5)
#' res3<-global.test(db_long, db_wide=NULL, from=1, to=6, quantiles=NULL, 
#'       times=times, tmat = tmat, replicas = 10, positions=positions,
#'       namesStates=namesStates, timesNames=timesNames, status=status)
#' 
#' library(mstate)   
#' data(prothr)
#' 
#' 
#' res4<-global.test(db_long = prothr, from=1, to=3, replicas=5, limit=0.90, 
#'         quantiles=c(.05, .10, .20, .30, 0.40))
#'         
#' round(res4$localTests,4)
#'
#' @author Gustavo Soutinho and Luis Meira-Machado.
#'
#' @importFrom "survidm" tprob markov.test 
#' @importFrom "survival" coxph Surv survfit strata untangle.specials
#' @importFrom "graphics" legend abline axis legend lines matplot par plot polygon
#' @importFrom "stats" pchisq pnorm quantile sd na.omit terms approxfun as.formula rnorm rpois weighted.mean
#' @importFrom "utils" capture.output
#' @importFrom "mstate" transMat msprep msfit probtrans LMAJ xsect cutLMms msdata2etm trans2Q events
#' @importFrom "stats" model.matrix model.frame model.response model.offset  
#' @importFrom "stats" delete.response delete.response

#' @export global.test
#' @export local.test
#' @export LR.tests
#' @export Cox.globaltest
#' @export prepMSM
#' @export transMatMSM
#' @export summary.markovMSM
#' @export print.markovMSM
#' @export plotMSM
#' @export eventsMSM
#' @S3method summary markovMSM
#' @S3method print markovMSM

global.test<- function(db_long, db_wide=NULL, from=1, to=3, quantiles=c(.05, 
                      .10, .20, .30, 0.40), times=NULL, 
                       tmat=NULL, replicas=10, limit=0.90,
                       positions=list(c(2, 3), c(3), c()),
                       namesStates =  c("Alive", "Rec",  "Death"),
                       timesNames = c(NA, "time1","Stime"),
                       status=c(NA, "event1","event")){
  
  
  i<-from
  
  #i<-2
  
  j<-to
  
  #j<-5
  
  M<-replicas 
  
  DIF<-NULL
  
  areaAll<-NULL
  
  ETAll<-NULL
  
  p.valueAll<-NULL
  
  ET_bootAll<-NULL
  
  ET.qiAll<-NULL
  
  ET.qsAll<-NULL
  
  ET.qi2All<-NULL
  
  ET.qs2All<-NULL
  
  times<-sort(times)
  
  if(is.null(quantiles) & is.null(times)){
    
    stop("One of the arguments quantiles or times must be NULL.")
  }
  
  if(!is.null(quantiles) & !is.null(times)){
    
    stop("Only one of the arguments quantiles or times must be NULL.")
  }
  
  if(length(times)!=5 & !is.null(times)){
    
    stop("The length of times must be 5.")
    
  }
  
  #min e quantile 0.75
  
  if(!is.null(times)){
    
    if(!is.null(db_wide)){
      
      if(length(namesStates)==3){
        
        timeToQuantiles<-timesNames[2]
        
        min<-min(as.numeric(db_wide[, timeToQuantiles]))
        
        val0.75<-as.numeric(quantile(db_wide[, timeToQuantiles], 0.75))
        
        
      }else{#inicio difere de 3
        
        if(i==1){#from
          
          #head(db_wide)
          
          valToQuantiles<-rep(0,nrow(db_wide))
          
          for(a1 in 1:nrow(db_wide)){
            
            #a1<-1
            
            valToQuantiles[a1]<-min(db_wide[a1, timesNames[-1]])
            
          }
          valToQuantiles<-valToQuantiles[valToQuantiles!=0] 
          
        }else{
          
          valToQuantiles<-db_wide[,timesNames[i]]
          valToQuantiles<-valToQuantiles[valToQuantiles!=0] 
        }
        
        min<-min(as.numeric(valToQuantiles))
        
        val0.75<-as.numeric(quantile(valToQuantiles, 0.75))
        
      }#fim de casos em que length(names) difere de 3
      
    }else{#caso db_wide=NULL
      
      tempQ<-unique(db_long$Tstop)
      
      min<-min(as.numeric(tempQ))
      
      val0.75<-as.numeric(quantile(tempQ, 0.75))
      
    }
    
    if(length(which(times<val0.75 & times>min & !is.null(times)))!=5){
      
      stop("All the times must be greater than the minumum time and the third quantile of time.")
      
    }
    
  }
  
  #
  
  if(is.null(times))
    
    times.to.test<-length(quantiles)
  
  if(is.null(quantiles))
    
    times.to.test<-length(times)
  
  for(h in 1:times.to.test){
    
    #h<-1
    
    print(h)
    
    #for(i3 in 1:length(colnames(db_wide))){
    
    #i3<-1
    
    #  if (colnames(db_wide)[i3]==timeToQuantiles)
    
    #    indexS<-i3
    #}
    #s<-as.numeric(quantile(db_wide[,indexS], quantiles))[h]
    
    
    if(!is.null(db_wide)){
      
      if(length(namesStates)==3){
        
        if(is.null(times)){
          
          timeToQuantiles<-timesNames[2]
          
          s<-as.numeric(quantile(db_wide[, timeToQuantiles], quantiles))[h]
          
          s_quant<-as.numeric(quantile(db_wide[, timeToQuantiles], quantiles))
          
        }else{
          
          
          s<-as.numeric(times)[h]
          s_quant<-as.numeric(times)
          
        }
        
      }else{ #inicio difere de 3
        
        if(is.null(times)){
          if(i==1){#from
            
            #head(db_wide)
            
            valToQuantiles<-rep(0,nrow(db_wide))
            
            for(a1 in 1:nrow(db_wide)){
              
              #a1<-1
              
              valToQuantiles[a1]<-min(db_wide[a1, timesNames[-1]])
              
            }
            valToQuantiles<-valToQuantiles[valToQuantiles!=0] 
            
          }else{
            
            valToQuantiles<-db_wide[,timesNames[i]]
            valToQuantiles<-valToQuantiles[valToQuantiles!=0] 
          }
          
          s<-as.numeric(quantile(valToQuantiles, quantiles))[h]
          
          s_quant<-as.numeric(quantile(valToQuantiles, quantiles))
          
        }else{
          s<-as.numeric(times)[h]
          s_quant<-as.numeric(times)
          
        }
        
        
      }#fim de casos em que length(names) difere de 3
      
    }else{#caso db_wide=NULL
      
      
      if(is.null(times)){
        
        #head(db_long)
        
        tempQ<-unique(db_long$Tstop)
        
        #summary(tempQ)
        
        s<-as.numeric(quantile(tempQ, quantiles))[h]
        
        tmat <- attr(db_long, "trans")
        
        s_quant<-as.numeric(quantile(tempQ, quantiles))
        
      }else{
        
        s<-as.numeric(times)[h]
        s_quant<-as.numeric(times)
        
      }
      
    }
    
    #h<-2
    
    #s<-c(10, 30, 100, 200, 300, 650)[h]  #10 erro
    
    c0 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = db_long)
    
    msf0 <- msfit(object = c0, vartype = "greenwood", trans = tmat)
    
    resAJ<-probtrans(msf0, predt=s)
    
    #resAJ[[i]]
    
    tempos <- resAJ[[1]][,1]
    t.limit <- quantile(tempos, limit)
    tempos <- tempos[tempos < t.limit]
    
    #trans<-as.numeric(na.omit(as.numeric(resAJ$trans)))
    
    #resAJ[[i]][1:10,]
    
    LMpt0 <- LMAJ2(db_long, s=s, from=i)
    
    #LMpt0[200:250,]
    
    #length(LMpt0$time)
    
    #tempos[tempos %in% LMpt0$time]
    
    #LMpt0$time[LMpt0$time %in% tempos]
    
    #tempos[tempos %in% LMpt0$time] == LMpt0$time[LMpt0$time %in% tempos]
    
    tComuns<-tempos[tempos %in% LMpt0$time]
    
    resAJf<-resAJ[[i]][resAJ[[i]][,1] %in% tComuns,]
    
    LMpt0f<-LMpt0[LMpt0$time %in% tComuns,]
    
    #head(resAJf[,2:(j+1)])
    
    #head(LMpt0f[,2:(j+1)])
    
    pr.ini<-resAJf[,2:(j+1)]-LMpt0f[,2:(j+1)]
    
    pr<-pr.ini[-1,] #; pr
    
    #head(pr)
    
    tempos2 <- diff(tComuns)
    
    #length(tempos2)
    
    #tempos2[1]*pr[1,]
    #tempos2[2]*pr[2,]
    
    #(pr * tempos2)[1:2,]
    
    area<-colSums((pr * tempos2))
    
    areaAll<-rbind(areaAll, cbind(s, as.data.frame(t(area))))
    
    #areaAll
    
    areasBoot<-matrix(NA, nrow = M, ncol =  j) #ncol(resAJ$trans)
    
    prTimes1<-cbind(tComuns, pr.ini)
    
    DIF<-rbind(DIF, cbind(s, prTimes1))
    
    #none<-NULL  (para filtro...)
    
    #for(n1 in from:to){
    #  #n1<-1
    #  if(n1==from)
    #    none<-c(none,unique(prTimes1[,(1+n1)]==resAJf[,(1+n1)]-1))
    #  else
    #    none<-c(none, unique(prTimes1[,(1+n1)]==resAJf[,(1+n1)])) 
    
    #}
    
    ET_boot<-NULL
    
    for(k in 1:M){
      
      #k<-1
      
      cat("Monte Carlo sample  =", k,"\n")
      
      #head(db_wide)
      
      #unique(db_wide$id)
      
      if(!is.null(db_wide)){
        
        n1 <- dim(db_wide)[1]
        pos <- sample.int(n1, n1, replace=TRUE)
        db_wide2 <- db_wide[pos,]
        
        
        #head(db_wide2)
        #positions<-list(c(2, 3), c(3), c())
        #namesStates = names = c("Alive", "Rec",  "Death")
        
        tmat <-transMatMSM(positions, namesStates)
        
        #times = c(NA, "time1","Stime")
        #status=c(NA, "event1","event")
        
        trans = tmat
        
        db_long2<- prepMSM(data=db_wide2, trans, timesNames, status)
        
        
        #db_long2<- msprep(data = db_wide2, trans = tmat, 
        #                 time = c(NA, "rec", "ae","recae", "rel", "srv"), 
        #                status = c(NA, "rec.s", "ae.s", "recae.s","rel.s", "srv.s"), 
        #                keep = c("match", "proph", "year", "agecl"))
        
        
      }else{
        
        
        n1<-length(unique(db_long$id))
        #n1 <- dim(db_wide)[1]
        pos <- sample.int(n1, n1, replace=TRUE)
        
        #table(pos)
        
        pos2<-sort(pos)
        
        db_long2<-NULL
        
        for(i4 in 1:n1){
          
          db_long2 <- rbind(db_long2,db_long[db_long$id %in% pos2[i4],])
          
        }
        
        #db_long2[db_long2$id==3,]
      }
      
      c0 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = db_long2)
      
      
      msf0 <- msfit(object = c0, vartype = "greenwood", trans = tmat)
      
      resAJ<-probtrans(msf0, predt=s)
      
      tempos <- resAJ[[1]][,1]
      
      t.limit <- quantile(tempos, limit)
      
      tempos <- tempos[tempos < t.limit]
      
      
      #
      
      msdata<-db_long2
      xss <- xsect(msdata, s)
      infrom <- xss$id[xss$state %in% from]
      msdatas <- cutLMms(msdata, LM = s)
      msdatasfrom <- msdatas[msdatas$id %in% infrom, ]
      
      #print(nrow(msdatasfrom))
      
      
      if(nrow(msdatasfrom)<2){
        
        
        print('Error')
        
      }else{
        
        
        LMpt0 <- LMAJ2(db_long2, s=s, from=i)
        
        #tempos[tempos %in% LMpt0$time] == LMpt0$time[LMpt0$time %in% tempos]
        
        tComuns<-tempos[tempos %in% LMpt0$time]
        
        resAJf<-resAJ[[i]][resAJ[[i]][,1] %in% tComuns,]
        
        LMpt0f<-LMpt0[LMpt0$time %in% tComuns,]
        
        pr.ini<-resAJf[,2:(j+1)]-LMpt0f[,2:(j+1)]
        
        pr<-pr.ini[-1,]; pr
        
        tempos2b <- diff(tComuns)
        
        areasBoot[k,]<-colSums((pr * tempos2b))
        
        #para o grafico
        
        prTimes<-cbind(tComuns, pr.ini)
        
        prTimes1$tComuns %in%  prTimes$tComuns
        
        
        dataframe<-as.data.frame(prTimes1$tComuns)
        
        colnames(dataframe)<-'tComuns'
        
        ET_boot_rep<-merge(x=dataframe, y=prTimes, by='tComuns',all.x = T)
        
        ET_boot<-rbind(ET_boot,cbind(ET_boot_rep,rep(k,nrow(ET_boot_rep)),  rep(s,nrow(ET_boot_rep))))
        
        #ET_boot[200:350,]
        
        #cat(k, '\n')
        
        #cat(h, '\n')
        
      }
      
      #
    } 
    
    sd_areasBoot<-rep(0, j) # ncol(resAJ$trans)
    
    for(l in 1: j){#ncol(resAJ$trans)
      
      #l<-1
      
      sd_areasBoot[l]<-sd(areasBoot[,l], na.rm = T)
      
    }
    
    
    ET<-area/sd_areasBoot; ET
    
    ETAll<-rbind(ETAll, cbind(s, as.data.frame(t(ET))))
    
    
    #none<-unique(none)
    
    #if(length(none)==1){
    
    #  p.value <- 2*(1-pnorm(abs(ET))); p.value
    
    #  for(n2 in from:to){
    
    #    p.value[n2]<-'None'
    
    #  }
    
    #}else{
    
    #  p.value <- 2*(1-pnorm(abs(ET))); p.value
    
    #}
    
    p.value <- 2*(1-pnorm(abs(ET))); p.value
    
    p.valueAll<-rbind(p.valueAll, cbind(s, as.data.frame(t(p.value ))))
    
    p.valueAll
    
    #
    
    colnames(ET_boot)[(j+2)]<-'replica'
    
    colnames(ET_boot)[(j+3)]<-'s'
    
    #head(ET_boot)
    
    len<-unique(ET_boot$tComuns)
    
    ET.qi<-matrix(NA, nrow = length(len), ncol = (j+1))
    
    ET.qs<-matrix(NA, nrow = length(len), ncol = (j+1))
    
    for(m in 1:length(len)){
      
      #m<-100
      
      #ET_boot[ET_boot$tComuns==len[m],]
      
      for(n in 2:(j+1)){
        #n<-2
        #ET_boot[ET_boot$tComuns==len[m],2]
        ET.qi[m,n]<-as.numeric(quantile(ET_boot[ET_boot$tComuns==len[m],n],probs =  0.025,na.rm = T))
        ET.qs[m,n]<-as.numeric(quantile(ET_boot[ET_boot$tComuns==len[m],n],probs = .975,na.rm = T))
        
      }
      
    }
    
    ET.qi[,1]<-len
    
    ET.qs[,1]<-len
    
    ET_boot[1:20,] 
    
    ET_boot2<-ET_boot
    
    for(p in 1:M){
      
      #p<-1
      
      #ET_boot[ET_boot$replica==p,]
      
      for(q in 2:(j+1)){
        
        #q<-2
        #ET_boot[ET_boot$replica==p,q]
        
        for (r in 1:length(len)){
          
          #r<-14
          ET_boot2[ET_boot2$replica==p,q][r]<-ifelse(is.na(ET_boot2[ET_boot2$replica==p,q][r]), 
                                                     ET_boot2[ET_boot2$replica==p,q][(r-1)],
                                                     ET_boot2[ET_boot2$replica==p,q][r])
          
          
        }
      }
    }
    
    #ET_boot2[ET_boot2$replica==p,]
    
    ET.qi2<-matrix(NA, nrow = length(len), ncol = (j+1))
    
    ET.qs2<-matrix(NA, nrow = length(len), ncol = (j+1))
    
    for(m2 in 1:length(len)){
      
      #m<-100
      
      #ET_boot[ET_boot$tComuns==len[m],]
      
      for(n2 in 2:(j+1)){
        #n<-2
        #ET_boot2[ET_boot2$tComuns==len[m],4]
        ET.qi2[m2,n2]<-as.numeric(quantile(ET_boot2[ET_boot2$tComuns==len[m2],n2],probs =  0.025,na.rm = T))
        ET.qs2[m2,n2]<-as.numeric(quantile(ET_boot2[ET_boot2$tComuns==len[m2],n2],probs = .975,na.rm = T))
        
      }
      
    }
    
    ET.qi2[,1]<-len
    
    ET.qs2[,1]<-len
    
    ET.qiAll<-rbind(ET.qiAll,cbind(as.data.frame(ET.qi), rep(s,length(len))))
    
    ET.qi2All<-rbind(ET.qi2All,cbind(as.data.frame(ET.qi2), rep(s,length(len))))
    
    ET.qsAll<-rbind(ET.qsAll,cbind(as.data.frame(ET.qs), rep(s,length(len))))
    
    ET.qs2All<-rbind(ET.qs2All,cbind(as.data.frame(ET.qs2), rep(s,length(len))))
    
  }#fim times
  
  p.value.f<-p.valueAll
  
  #round(p.value.f,4)
  
  
  minAll<-rep(0, j)
  
  if(is.null(times)){
    
    minInt<-rep(0, (length(quantiles)-1))
    
  }else{
    
    minInt<-rep(0, (length(times)-1))
  }
  
  
  for (i2 in 1:j){
    
    #i2<-1
    
    for(j2 in 1: (nrow(p.valueAll)-1)){
      
      #j2<-1
      
      p.valueAll[,(i2+1)]<-ifelse(p.valueAll[,(i2+1)]== 'NaN', NA, p.valueAll[,(i2+1)])
      
      minInt[j2]<-mean(p.valueAll[j2:(j2+1),(i2+1)],na.rm = T)
      
    }
    
    #is.na(minInt)
    
    minAll[i2]<-min(minInt, na.rm=T)
  }
  
  minAll
  
  #head(ET.qiAll)
  
  names(ET.qiAll)[length(names(ET.qiAll))]<-'s'
  
  names(ET.qiAll)[1]<-'tComuns'
  
  names(ET.qiAll)[2:(length(3:ncol(DIF))+1)]<-names(DIF)[3:ncol(DIF)]
  
  #
  
  #head(ET.qi2All)
  
  names(ET.qi2All)[length(names(ET.qi2All))]<-'s'
  
  names(ET.qi2All)[1]<-'tComuns'
  
  names(ET.qi2All)[2:(length(3:ncol(DIF))+1)]<-names(DIF)[3:ncol(DIF)]
  
  #
  
  #head(ET.qsAll)
  
  names(ET.qsAll)[length(names(ET.qsAll))]<-'s'
  
  names(ET.qsAll)[1]<-'tComuns'
  
  names(ET.qsAll)[2:(length(3:ncol(DIF))+1)]<-names(DIF)[3:ncol(DIF)]
  
  #
  
  #head(ET.qs2All)
  
  names(ET.qs2All)[length(names(ET.qs2All))]<-'s'
  
  names(ET.qs2All)[1]<-'tComuns'
  
  names(ET.qs2All)[2:(length(3:ncol(DIF))+1)]<-names(DIF)[3:ncol(DIF)]
  
  
  if(!is.null(db_wide)){
    
    p.value.f2<-p.value.f[,c(1, (from+1):ncol(p.value.f))]
    
    minAllf<-minAll[from:length(minAll)]
    
  }else{
    
    p.value.f2<-p.value.f
    
    minAllf<-minAll
  }
  
  if(is.null(times)){
    res <-list(globalTest=minAllf, localTests=p.value.f2, 
               quantiles=as.numeric(s_quant), times=NULL, DIF=DIF, from=from, to=to,
               ET.qiAll=ET.qiAll, ET.qi2All=ET.qi2All, ET.qsAll=ET.qsAll, ET.qs2All=ET.qs2All,
               replicas=M, limit=limit)
  }else{
    
    res <-list(globalTest=minAllf, localTests=p.value.f2, 
               quantiles=NULL,times=times, DIF=DIF, from=from, to=to,
               ET.qiAll=ET.qiAll, ET.qi2All=ET.qi2All, ET.qsAll=ET.qsAll, ET.qs2All=ET.qs2All,
               replicas=M, limit=limit)
  }
  
  
  
  class(res) <- c("globalTest", "markovMSM")
  
  res$call <- match.call()
  
  return(res)
  
}

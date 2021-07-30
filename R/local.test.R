#' Area Under the two curves (AUC) local test for checking the Markov condition on Multi-state Models.
#'
#' @description This function is used to obtain the local test based on the AUC Markov Test for selected 
#' times and transitions between states.
#'
#' @param db_long A dataframe in long format containing the subject \code{id}; \code{from} corresponding 
#' to the starting state; the receiving state, \code{to}; the transition number, \code{trans}; the starting
#' time of the transition given by \code{Tstart}; the stopping time of the transition, \code{Tstop}, and
#' \code{status} for the  status variable, with 1 indicating an event (transition), 0 a censoring. 
#' @param db_wide Data frame in wide format in which to interpret time, status, id or keep, if appropriate.
#' @param times The starting times of the transition probabilities for checking the local Markov assumption
#' in Multi-state models.
#' @param from The starting state of the transition probabilities.
#' @param to The last receiving state considered for the estimation of the transition probabilities. All 
#' the probabilities among the first and the last states are also computed.
#' @param tmat The transition matrix for multi-state model.
#' @param replicas Number of replicas for the Monte Carlo simulation to standardization of the T-statistic
#' given by the difference of the areas of AJ and LM transition probabilities estimates.
#' @param limit Percentile of the event time used as the upper bound for the computation of the AUC-based test. 
#' @param positions List of possible transitions; x[[i]] consists of a vector of state numbers reachable from state i.
#' @param namesStates A character vector containing the names of either the competing risks or the states in the 
#' multi-state model specified by the competing risks or illness-death model. names should have the same 
#' length as the list x (for transMat), or either K or K+1 (for trans.comprisk), or 3 (for trans.illdeath). 
#' @param timesNames Either 1) a matrix or data frame of dimension n x S (n being the number of individuals 
#' and S the number of states in the multi-state model), containing the times at which the states are 
#' visited or last follow-up time, or 2) a character vector of length S containing the column names 
#' indicating these times. In the latter cases, some elements of time may be NA, see Details
#' @param status Either 1) a matrix or data frame of dimension n x S, containing, for each of the states, 
#' event indicators taking the value 1 if the state is visited or 0 if it is not (censored), or 2) a character
#' vector of length S containing the column names indicating these status variables. In the latter cases, some
#' elements of status may be NA, see Details. 
#'
#' @return An object with a list with the following outcomes:
#' \item{localTest}{p-value of AUC local tests for each times and transitions.}
#' \item{trans}{The transition matrix describing the states and transitions of the multi-state model.}
#' \item{times}{Times selected for the AUC Local tests.}
#' \item{DIF}{Differences between the AJ and the LMAJ estimates for each transition probabilites from the 
#' starting state until the receiving state given by only one replica where 's' represent each of the 
#' quantile times.}  
#' \item{from}{The starting state considered for the AUC Local tests.}
#' \item{to}{The last receiving state considered for the the AUC Local tests.}
#' \item{ET.qiAll}{The lower limit of the diferences between the AJ and the LMAJ estimates given by the
#' Monte Carlo simulation in each transition for each "s" quantile times. \code{ET.qi2All} means the same
#' but missing values were replaced by the previous diferences of estimators.} 
#' \item{ET.qsAll}{The upper limit of the diferences between the AJ and the LMAJ estimates given by the Monte
#' Carlo simulation in each transition for each "s" quantile times. \code{ET.qs2All} means the same but
#' missing values were replaced by the previous diferences of estimators.} 
#' \item{replicas}{Number of replicas for the Monte Carlo simulation.}
#' \item{limit}{Percentil of the times used in the AUC local tests.}
#' \item{x}{x.}
#' \item{y}{y.}
#' 
#' 
#'
#' @examples
#' 
#' data("colonIDM")
#' colonIDM$event1[ colonIDM$time1==colonIDM$Stime & colonIDM$event1==1]<-0
#' db_wide<-colonIDM
#' 
#' positions<-list(c(2, 3), c(3), c())
#' namesStates =  c("Alive", "Rec",  "Death")
#' tmat <-transMatMSM(positions, namesStates)
#' timesNames = c(NA, "time1","Stime")
#' status=c(NA, "event1","event")
#' trans = tmat
#' db_long<- prepMSM(data=db_wide, trans, timesNames, status)
#' times=c(90, 180, 365, 730, 1095, 1460)
#' res<-local.test(db_long, db_wide, times=times, from=2, to=3, replicas=5, tmat = tmat)
#' res$localTest
#' 
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
#' times=c(90, 180, 365, 730, 1095, 1460)
#' res2<-local.test(db_long, db_wide, times=times, from=3, to=6, tmat = tmat, 
#'                replicas = 10, positions=positions, namesStates=namesStates, 
#'                timesNames=timesNames,status=status)
#' res2$localTest
#' 
#' 
#' library(mstate)
#' data(prothr)
#' 
#' res3<-local.test(db_long=prothr, times=times, from=2, to=3, replicas = 5, 
#' positions=positions, namesStates=namesStates, timesNames=timesNames,status=status)
#' round(res3$localTest, 4)
#'                
#' 
#' @author Gustavo Soutinho and Luis Meira-Machado.
#' 


local.test<- function(db_long, db_wide=NULL, times, from=1, to=3, tmat = NULL,
                      replicas=10, limit=0.90,
                      positions=list(c(2, 3), c(3), c()),
                      namesStates =  c("Alive", "Rec",  "Death"),
                      timesNames = c(NA, "time1","Stime"),
                      status=c(NA, "event1","event")){
  
  i<-from
  
  #i<-1
  
  j<-to
  
  #j<-6
  
  M<-replicas 
  
  DIF<-NULL
  
  areaAll<-NULL
  
  ETAll<-NULL
  
  p.valueAll<-NULL
  
  AJall<-NULL
  
  LMAJall<-NULL
  
  ET.qiAll<-NULL
  
  ET.qsAll<-NULL
  
  ET.qi2All<-NULL
  
  ET.qs2All<-NULL
  
  
  if(is.null(db_wide))
    tmat <- attr(db_long, "trans")
  
  
  for(h in 1:length(times)){
    
    #h<-2
    
    print(h)
    
    s<-times[h]
    
    c0 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = db_long)
    
    msf0 <- msfit(object = c0, vartype = "greenwood", trans = tmat)
    
    resAJ<-probtrans(msf0, predt=s)
    
    tempos <- resAJ[[1]][,1]
    
    t.limit <- quantile(tempos, 0.90)
    
    tempos <- tempos[tempos < t.limit]
    
    #trans<-as.numeric(na.omit(as.numeric(resAJ$trans)))
    
    resAJ[[1]] #from 1 to 6
    
    #resAJ[[i]][1:10,1:8]
    
    LMpt0 <- LMAJ2(db_long, s=s, from=i)
    
    #LMpt0[1:40,1:8]
    
    length(LMpt0$time)
    
    tempos[tempos %in% LMpt0$time]
    
    LMpt0$time[LMpt0$time %in% tempos]
    
    tempos[tempos %in% LMpt0$time] == LMpt0$time[LMpt0$time %in% tempos]
    
    tComuns<-tempos[tempos %in% LMpt0$time]
    
    resAJf<-resAJ[[i]][resAJ[[i]][,1] %in% tComuns,]
    
    LMpt0f<-LMpt0[LMpt0$time %in% tComuns,]
    
    pr.ini<-resAJf[,2:(j+1)]-LMpt0f[,2:(j+1)]
    
    pr<-pr.ini[-1,]#; pr
    
    #head(pr)
    
    pr<-pr[-1,]
    
    #head(pr)
    
    tempos2 <- diff(tComuns)
    
    #length(tempos2)
    
    #tempos2[1]*pr[1,]
    #tempos2[2]*pr[2,]
    
    #pr[1:10,]
    #tempos2[1:10]
    #(pr * tempos2)[1:10,]
    
    area<-colSums((pr * tempos2), na.rm = T) 
    
    areaAll<-rbind(areaAll, cbind(s, as.data.frame(t(area))))
    
    #areaAll
    
    areasBoot<-matrix(NA, nrow = M, ncol =  j) #ncol(resAJ$trans)
    
    prTimes1<-cbind(tComuns, pr.ini)
    
    DIF<-rbind(DIF, cbind(s, prTimes1))
    
    AJall<-rbind(AJall, cbind(s, resAJf[1:(j+1)]))
    
    LMAJall<-rbind(LMAJall, cbind(s, LMpt0f[1:(j+1)]))
    
    
    none<-NULL
    
    for(n1 in from:to){
      #n1<-1
      if(n1==from)
        none<-c(none,unique(prTimes1[,(1+n1)]==resAJf[,(1+n1)]-1))
      else
        none<-c(none, unique(prTimes1[,(1+n1)]==resAJf[,(1+n1)])) 
      
    }
    
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
      t.limit <- quantile(tempos, 0.90)
      tempos <- tempos[tempos < t.limit]
      
      LMpt0 <- LMAJ2(db_long2, s=s, from=i)
      
      #tempos[tempos %in% LMpt0$time] == LMpt0$time[LMpt0$time %in% tempos]
      
      tComuns<-tempos[tempos %in% LMpt0$time]
      
      resAJf<-resAJ[[i]][resAJ[[i]][,1] %in% tComuns,]
      
      #resAJf<-resAJ[[1]][resAJ[[1]][,1] %in% tComuns,] 
      
      LMpt0f<-LMpt0[LMpt0$time %in% tComuns,]
      
      #o que estava 
      
      #pr<-resAJf[,2:(j+1)]-LMpt0f[,2:(j+1)]
      
      #pr<-pr[-1,]; pr
      
      #tempos2 <- diff(tComuns)
      
      #areasBoot[k,]<-colSums((pr * tempos2), na.rm = T) #,
      
      
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
      
      
    } 
    
    sd_areasBoot<-rep(0, j) # ncol(resAJ$trans)
    
    
    for(l in 1: j){#ncol(resAJ$trans)
      
      #l<-3
      
      sd_areasBoot[l]<-sd(areasBoot[,l], na.rm = T)
      
    }
    
    
    ET<-area/sd_areasBoot; ET
    
    ETAll<-rbind(ETAll, cbind(s, as.data.frame(t(ET))))
    
    #ETAll
    
    
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
    
    ET_boot[1:20,] #with missing values. how to replace with the previous value
    
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
    
  }#fim tempos
  
  
  if(!is.null(db_wide)){
    
    p.valueAll.f<-p.valueAll[,c(1, (from+1):ncol(p.valueAll))]
    
  }else{
    
    
    p.valueAll.f<-p.valueAll
    
  }
  round(p.valueAll.f,4)
  
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
  
  
  res<-list(localTest=p.valueAll.f, trans = tmat, times=times, DIF=DIF, AJall=AJall, LMAJall=LMAJall, from=from, 
            to=to, ET.qiAll=ET.qiAll, ET.qi2All=ET.qi2All, ET.qsAll=ET.qsAll, ET.qs2All=ET.qs2All, replicas=M, 
            limit=limit)
  
  
  class(res) <- c("localTest", "markovMSM")
  
  res$call <- match.call()
  
  return(res)
  
}
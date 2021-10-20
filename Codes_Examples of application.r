library(markovMSM)
data("colonMSM")
db_wide<-colonMSM
#db_wide[db_wide$event1==1 & db_wide$event==1 & db_wide$time1==db_wide$Stime,2]<-0
head(db_wide[,1:11])
positions<-list(c(2, 3), c(3), c())
namesStates = c("Alive", "Rec", "Death")
tmat <-transMatMSM(positions, namesStates)
timesNames = c(NA, "time1","Stime")
status=c(NA, "event1","event")
trans = tmat
db_long<- prepMSM(data=db_wide, trans, timesNames, status)
db_long[1:10,]
eventsMSM(db_long)

#library(survidm)
#nevents(colonIDM)
res<-Cox.globaltest(data=db_long, transition = 3)
res

summary(res$p.value)

print(res)

res2<-local.test(db_long, db_wide, times=180, from=1, to=3, replicas=100,
                 tmat = tmat)
res2$localTest

res3<-local.test(db_long, db_wide, times=180, from=2, to=3, replicas=100,
                 tmat = tmat)
res3$localTest

plot(res2, tran='pstate2', axis.scale=c(0,0.25), difP=FALSE)
plot(res3, tran='pstate3', axis.scale=c(0,1), difP=FALSE)

plot(res2, tran='pstate2',axis.scale=c(-0.03,0.03), difP=TRUE)
plot(res3, tran='pstate3', axis.scale=c(-0.30,.10), difP=TRUE)

res4<-global.test(db_long, db_wide, from=1, to=3, replicas = 10, tmat=tmat)
round(res4$globalTest,3)

res5<-global.test(db_long, db_wide, from=2, to=3, replicas = 10, tmat=tmat)
round(res5$globalTest,3)

round(res4$localTest,3)

round(res5$localTest,3)

res5b<-global.test(db_long, db_wide, from=2, to=3, replicas = 10, 
                   quantiles=c(0.05, .1, .20, .30, .40),times=c(102.4, 173.0, 290.6, 469.2, 726.8),
                   tmat=tmat) #confirma-se erro

res5c<-global.test(db_long, db_wide, from=2, to=3, replicas = 10, 
                   times=c(102.4, 173.0, 290.6, 469.2, 726.8),
                   tmat=tmat) #confirma-se erro

res5c<-global.test(db_long, db_wide, from=2, to=3, replicas = 10, quantiles=NULL,
                   times=c(102.4, 173.0, 290.6, 469.2, 726.8),
                   tmat=tmat)

round(res5c$globalTest,3) 

round(res5c$localTest,3)


res5d<-global.test(db_long, db_wide, from=2, to=3, replicas = 10, quantiles=NULL,times=NULL, tmat=tmat)

plot(res4, quantileOrder=3, axis.scale=c(-.04, .02))

plot(res5, quantileOrder=3, axis.scale=c(-.10, .20))

res6<-LR.tests(db_long=db_long, times=180, transition = 3, replicas = 5)
res6$localTestLR

times<-c(102.4, 173, 290.6, 469.2, 726.8)
res7<-LR.tests(db_long=db_long, times=times, transition = 3, replicas = 5)
res7$localTestLR

res7$globalTestLR


data("ebmt4")
db_wide <- ebmt4
positions=list(c(2, 3, 5, 6), c(4, 5, 6), c(4, 5, 6),
               c(5, 6), c(6), c())
namesStates = c("Tx", "Rec", "AE", "Rec+AE", "Rel", "Death")
tmat <-transMatMSM(positions, namesStates)
timesNames = c(NA, "rec", "ae","recae", "rel", "srv")
status=c(NA, "rec.s", "ae.s", "recae.s","rel.s", "srv.s")
trans = tmat
db_long<- prepMSM(data=db_wide, trans, timesNames, status)
db_long[1:10,]

res8<-global.test(db_long, db_wide, from=1, to=5,
                  quantiles=c(.05, .10, .20, .30, 0.40), times=c(9.9, 12, 15, 18, 21),
                  tmat = tmat, replicas = 10,
                  positions=positions, namesStates=namesStates,
                  timesNames=timesNames, status=status)

res8<-global.test(db_long, db_wide, from=1, to=5, times=c(9.9, 12, 15, 18, 21),
                  tmat = tmat, replicas = 100,
                  positions=positions, namesStates=namesStates,
                  timesNames=timesNames, status=status)


res8<-global.test(db_long, db_wide, from=1, to=5,
                  quantiles=c(.05, .10, .20, .30, 0.40),
                  tmat = tmat, replicas = 10,
                  positions=positions, namesStates=namesStates,
                  timesNames=timesNames, status=status)

round(res8$globalTest,4)
round(res8$localTests,4)

res9<-global.test(db_long = prothr, db_wide = NULL, from=2, to=3,
                  replicas=10, limit=0.90,
                  quantiles=c(.05, .10, .20, .30, 0.40))


round(res9$globalTest,4)


times<-c(73.5, 117 ,223, 392, 681)
res10<-LR.tests(db_long=prothr, times=times, transition = 4, replicas = 10)
res10$localTestLR

res10$globalTestLR

res11<-LR.tests(db_long=prothr, times=times, transition = 3, replicas = 10)
res11$localTestLR

res11$globalTestLR




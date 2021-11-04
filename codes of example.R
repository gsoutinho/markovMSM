library(markovMSM)
data("colonMSM")
db_wide <- colonMSM
head(db_wide[c(1:2,16,21),1:11])
positions <-list(c(2, 3), c(3), c())
state.names <- c("Alive", "Rec", "Death")
tmat <-transMatMSM(positions, state.names)
timesNames <- c(NA, "time1", "Stime")
status <- c(NA, "event1", "event")
trans <- tmat
db_long <- prepMSM(data = db_wide, trans, timesNames, status)
db_long[1:10,]
eventsMSM(db_long)

res <- PHM.test(data=db_long, from=2, to=3)
res

set.seed(1234)
res2<-AUC.test(db_long, db_wide, times=180, from=1, to=3, type='local',
               replicas=100, tmat = tmat)

res2$localTest

set.seed(1234)
res3<-AUC.test(db_long, db_wide, times=180, from=2, to=3, type='local',
               replicas=100, tmat = tmat)
res3$localTest

plot(res2, to=2, axis.scale=c(0,0.25), difP=FALSE)
plot(res3, to=3, axis.scale=c(0,1), difP=FALSE)

plot(res2, to=2, axis.scale=c(-0.03,0.03), difP=TRUE)
plot(res3, to=3, axis.scale=c(-0.30,.10), difP=TRUE)


set.seed(1234)
res4<-AUC.test(db_long, db_wide, from=1, to=3, replicas = 100, tmat=tmat)
round(res4$globalTest,3)

set.seed(1234)
res5<-AUC.test(db_long, db_wide, from=2, to=3, type='global', replicas = 100,
               tmat=tmat)
round(res5$globalTest,3)

round(res4$localTest,3)

round(res5$localTest,3)

plot(res4, quantileOrder=2, axis.scale=c(-.04, .02))
plot(res5, quantileOrder=2, axis.scale=c(-.10, .20))

set.seed(1234)

res6<-LR.test(db_long=db_long, times=180, from = 2, to = 3, replicas = 1000)
res6$globalTestLR


data("ebmt4")
db_wide <- ebmt4
positions=list(c(2, 3, 5, 6), c(4, 5, 6), c(4, 5, 6),
               c(5, 6), c(6), c())
namesStates = c("Tx", "Rec", "AE", "Rec+AE", "Rel", "Death")
tmat <-transMatMSM(positions, namesStates)
timesNames = c(NA, "rec", "ae","recae", "rel", "srv")
status=c(NA, "rec.s", "ae.s", "recae.s","rel.s", "srv.s")
trans = tmat
db_long <- prepMSM(data=db_wide, trans, timesNames, status)
db_long[1:10,]

set.seed(1234)
res7<-AUC.test(db_long, db_wide, from=1, to=5, type='global',
               quantiles=c(.05, .10, .20, .30, 0.40),
               tmat = tmat, replicas = 100,
               positions=positions, namesStates=namesStates,
               timesNames=timesNames, status=status)

round(res7$globalTest, 4)

round(res7$localTests,4)

set.seed(1234)
res8<-AUC.test(db_long = prothr, db_wide = NULL, from=2, to=3,
               type='global', replicas=100, limit=0.90,
               quantiles=c(.05, .10, .20, .30, 0.40))

round(res8$globalTest,5)

round(res8$localTests,4)

set.seed(1234)
times<-c(73.5, 117, 223, 392, 681)
res9<-LR.test(db_long=prothr, times=times, from = 2, to = 3, replicas = 1000)

res9$localTestLR
res9$globalTestLR

set.seed(1234)
res10<-LR.test(db_long=prothr, times=times, from = 2, to = 1, replicas = 1000)
res10$localTestLR
res10$globalTestLR

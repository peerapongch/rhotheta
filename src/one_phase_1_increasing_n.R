library(grid)
library(gridExtra)
library(caret)
library(mboost)
library(plyr)
library(ggplot2)
library(reshape2)
source('../functions_plot.R')
setwd("C:/cs350/models/sprint1")
mse <- function(ytrue,yhat){
  return(mean((ytrue-yhat)^2))
}


## n=100
nsam <- 100
nrep <- 100
nuisance<-'pure'
load(paste('../../data/rdas/traindata_',nsam,'_',nrep,'_',nuisance,'.rda',sep=''))

m_gam1 <- gamboost(theta~.,
                   data=train,
                   control=boost_control(mstop=500,
                                         trace=TRUE
                   ),
                   baselearner = "bbs",
                   family=Gaussian(),
                   dfbase=4
)
system.time(
  cvm1 <- cvrisk(m_gam1, folds = cv(weights = model.weights(m_gam1),
                                    type = "kfold"))
)
m_gam1[mstop(cvm1)]

## n=200
nsam <- 200
nrep <- 100
nuisance<-'pure'
load(paste('../../data/rdas/traindata_',nsam,'_',nrep,'_',nuisance,'.rda',sep=''))

m_gam2 <- gamboost(theta~.,
                   data=train,
                   control=boost_control(mstop=500,
                                         trace=TRUE
                   ),
                   baselearner = "bbs",
                   family=Gaussian(),
                   dfbase=4
)
system.time(
  cvm2 <- cvrisk(m_gam2, folds = cv(weights = model.weights(m_gam2),
                                    type = "kfold"))
)
m_gam2[mstop(cvm2)]

## n=500
nsam <- 500
nrep <- 100
nuisance<-'pure'
load(paste('../../data/rdas/traindata_',nsam,'_',nrep,'_',nuisance,'.rda',sep=''))

m_gam3 <- gamboost(theta~.,
                   data=train,
                   control=boost_control(mstop=500,
                                         trace=TRUE
                   ),
                   baselearner = "bbs",
                   family=Gaussian(),
                   dfbase=4
)
system.time(
  cvm3 <- cvrisk(m_gam3, folds = cv(weights = model.weights(m_gam3),
                                    type = "kfold"))
)
m_gam3[mstop(cvm3)]


## n=1000
nsam <- 1000
nrep <- 100
nuisance<-'pure'
load(paste('../../data/rdas/traindata_',nsam,'_',nrep,'_',nuisance,'.rda',sep=''))

m_gam4 <- gamboost(theta~.,
                   data=train,
                   control=boost_control(mstop=500,
                                         trace=TRUE
                   ),
                   baselearner = "bbs",
                   family=Gaussian(),
                   dfbase=4
)
system.time(
  cvm4 <- cvrisk(m_gam4, folds = cv(weights = model.weights(m_gam4),
                                    type = "kfold"))
)
m_gam4[mstop(cvm4)]



## prediction
nsam <- 100
nrep <- 100
nuisance<-'pure'
load(paste('../../data/rdas/traindata_',nsam,'_',nrep,'_',nuisance,'.rda',sep=''))

pred_gam1<-predict(m_gam1,newdata=test.X)
mse1 <- mse(pred_gam1,test.Y)
bias1 <- mean(pred_gam1-test.Y)
var1 <- var(pred_gam1)
pred_gam1<-summarisePred(pred_gam1,test.Y)
p1<-plotPred(pred_gam1,'(a) Boosted GAM, n=100',20,15)


nsam <- 200
nrep <- 100
nuisance<-'pure'
load(paste('../../data/rdas/traindata_',nsam,'_',nrep,'_',nuisance,'.rda',sep=''))

pred_gam2<-predict(m_gam2,newdata=test.X)
mse2 <- mse(pred_gam2,test.Y)
bias2 <- mean(pred_gam2-test.Y)
var2 <- var(pred_gam2)
pred_gam2<-summarisePred(pred_gam2,test.Y)
p2<-plotPred(pred_gam2,'(b) Boosted GAM, n=200',20,15)

nsam <- 500
nrep <- 100
nuisance<-'pure'
load(paste('../../data/rdas/traindata_',nsam,'_',nrep,'_',nuisance,'.rda',sep=''))

pred_gam3<-predict(m_gam3,newdata=test.X)
mse3 <- mse(pred_gam3,test.Y)
pred_gam3<-summarisePred(pred_gam3,test.Y)
p3<-plotPred(pred_gam3,'(c) Boosted GAM, n=500',20,15)

nsam <- 1000
nrep <- 100
nuisance<-'pure'
load(paste('../../data/rdas/traindata_',nsam,'_',nrep,'_',nuisance,'.rda',sep=''))

pred_gam4<-predict(m_gam4,newdata=test.X)
mse4 <- mse(pred_gam4,test.Y)
pred_gam4<-summarisePred(pred_gam4,test.Y)
p4<-plotPred(pred_gam4,'(d) Boosted GAM, n=1,000',20,15)

gam_plots <- grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2)

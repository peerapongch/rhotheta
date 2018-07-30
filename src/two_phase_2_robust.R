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
scaleback <- function(data,mus,sigmas){
  df<-scale(data,center=FALSE,scale=1/sigmas)
  df<-scale(df,center=-mus,scale=FALSE)
  return(df)
}

## train boosting at 197 
nsam <- 197
nrep <- 100
nuisance<-'pure'
load(paste('../../data/rdas/traindata_',nsam,'_',nrep,'_',nuisance,'.rda',sep=''))
the_mus<-mus_train
the_sigmas<-sigmas_train

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

##test data: rho 40 
nsam <- 197
nrep <- 100
nuisance<-'rho40'
load(paste('../../data/rdas/traindata_',nsam,'_',nrep,'_',nuisance,'.rda',sep=''))

h <- sum(seq(1,nsam-1,1)^(-1))
pred_w <- (test.X$S_n*sigmas_test['S_n']+mus_test['S_n'])/h
pred_d <- test.X$Pi_n*sigmas_test['Pi_n']+mus_test['Pi_n']
eta1<-test.X$X1*sigmas_test['X1']+mus_test['X1']
mult<-(nsam-1)/nsam
pred_fu<-eta1*mult


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

pred_gam1<-predict(m_gam1,newdata=test.X)
pred_gam1<-summarisePred(pred_gam1,test.Y)
pred_w1 <- summarisePred(pred_w,test.Y)
pred_d1 <- summarisePred(pred_d,test.Y)
pred_fu1 <- summarisePred(pred_fu,test.Y)

real_test<-scaleback(test.X,mus_train,sigmas_train)
real_test<-data.frame(scale(real_test,center=the_mus,scale=the_sigmas))

pred_gam1_2<-predict(m_gam1,newdata=real_test)
pred_gam1_2<-summarisePred(pred_gam1_2,test.Y)
p1_2<-plotPred(pred_gam1_2,'test',20,15)

## test data: rho 80 
nsam <- 197
nrep <- 100
nuisance<-'rho80'
load(paste('../../data/rdas/traindata_',nsam,'_',nrep,'_',nuisance,'.rda',sep=''))

h <- sum(seq(1,nsam-1,1)^(-1))
pred_w <- (test.X$S_n*sigmas_test['S_n']+mus_test['S_n'])/h
pred_d <- test.X$Pi_n*sigmas_test['Pi_n']+mus_test['Pi_n']
eta1<-test.X$X1*sigmas_test['X1']+mus_test['X1']
mult<-(nsam-1)/nsam
pred_fu<-eta1*mult

pred_gam2<-predict(m_gam1,newdata=test.X)
pred_gam2<-summarisePred(pred_gam2,test.Y)
pred_w2 <- summarisePred(pred_w,test.Y)
pred_d2 <- summarisePred(pred_d,test.Y)
pred_fu2 <- summarisePred(pred_fu,test.Y)

real_test2<-scaleback(test.X,mus_train,sigmas_train)
real_test2<-data.frame(scale(real_test2,center=the_mus,scale=the_sigmas))

pred_gam2_2<-predict(m_gam2,newdata=real_test)
pred_gam2_2<-summarisePred(pred_gam2_2,test.Y)
p2_2<-plotPred(pred_gam2_2,'test',20,15)

## test data: rho 120
nsam <- 197
nrep <- 100
nuisance<-'rho120'
load(paste('../../data/rdas/traindata_',nsam,'_',nrep,'_',nuisance,'.rda',sep=''))

h <- sum(seq(1,nsam-1,1)^(-1))
pred_w <- (test.X$S_n*sigmas_test['S_n']+mus_test['S_n'])/h
pred_d <- test.X$Pi_n*sigmas_test['Pi_n']+mus_test['Pi_n']
eta1<-test.X$X1*sigmas_test['X1']+mus_test['X1']
mult<-(nsam-1)/nsam
pred_fu<-eta1*mult

pred_gam3<-predict(m_gam1,newdata=test.X)
pred_gam3<-summarisePred(pred_gam3,test.Y)
pred_w3 <- summarisePred(pred_w,test.Y)
pred_d3 <- summarisePred(pred_d,test.Y)
pred_fu3 <- summarisePred(pred_fu,test.Y)

real_test3<-scaleback(test.X,mus_train,sigmas_train)
real_test3<-data.frame(scale(real_test3,center=the_mus,scale=the_sigmas))

pred_gam3_2<-predict(m_gam3,newdata=real_test)
pred_gam3_2<-summarisePred(pred_gam3_2,test.Y)
p3_2<-plotPred(pred_gam3_2,'test',20,15)

### plots 
p1<-plotPred(pred_gam1,expression(paste('(a) Boosted GAM, ',rho,'=40',sep='')),15,10)
p2<-plotPred(pred_w1,expression(paste('(b) Watterson\'s Estimator, ',rho,'=40',sep='')),15,10)
p3<-plotPred(pred_d1,expression(paste('(c) Tajima\'s Estimator, ',rho,'=40',sep='')),15,10)
p4<-plotPred(pred_fu1,expression(paste('(d) Fu and Li\'s Estimator, ',rho,'=40',sep='')),15,10)

p5<-plotPred(pred_gam2,expression(paste('(e) Boosted GAM, ',rho,'=80',sep='')),15,10)
p6<-plotPred(pred_w2,expression(paste('(f) Watterson\'s Estimator, ',rho,'=80',sep='')),15,10)
p7<-plotPred(pred_d2,expression(paste('(g) Tajima\'s Estimator, ',rho,'=80',sep='')),15,10)
p8<-plotPred(pred_fu2,expression(paste('(h) Fu and Li\'s Estimator, ',rho,'=80',sep='')),15,10)

p9<-plotPred(pred_gam1,expression(paste('(i) Boosted GAM, ',rho,'=120',sep='')),15,10)
p10<-plotPred(pred_w3,expression(paste('(j) Watterson\'s Estimator, ',rho,'=120',sep='')),15,10)
p11<-plotPred(pred_d3,expression(paste('(k) Tajima\'s Estimator, ',rho,'=120',sep='')),15,10)
p12<-plotPred(pred_fu3,expression(paste('(l) Fu and Li\'s Estimator, ',rho,'=120',sep='')),15,10)

#plots<-grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,nrow=3,ncol=4)
plots<-grid.arrange(p1,p5,p9,p2,p6,p10,p3,p7,p11,p4,p8,p12,nrow=4,ncol=3)

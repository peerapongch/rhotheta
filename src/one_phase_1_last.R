library(caret)
library(mboost)
library(plyr)
library(ggplot2)
source('../functions_plot.R')
# =================== Section 0: Initialising ===================
setwd("C:/cs350/models/sprint1")
# change the input rda files to train different models and make different predictions
# this allows investigation into making prediction of datasets arising from different relationships of n and theta
nsam <- 100
nrep <- 70
load(paste('../../data/rdas/traindata_',nsam,'_',nrep,'_pure.rda',sep=''))

# =================== Section 1: default implementation ===================
# this setting uses penalised (k=2) cubic spline with 20 knots, df=4, and no boundary knots
# hence it lacks the linear extrapolation capabiltiy of a natural cublic smoothing splines
# fitting 
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
  cvm1 <- cvrisk(m_gam2, folds = cv(weights = model.weights(m_gam2),
                                    type = "kfold"))
)
m_gam2[mstop(cvm1)]

nsam <- 100
nrep <- 30
load(paste('../../data/rdas/traindata_',nsam,'_',nrep,'_pure.rda',sep=''))

pred_gam<-predict(m_gam2,newdata=test.X)
pred_gam<-summarisePred(pred_gam,test.Y)

p<-plotPred(pred_gam,'sometitle',20,15)

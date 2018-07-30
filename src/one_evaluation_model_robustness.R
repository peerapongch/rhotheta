library(plyr)
library(ggplot2)
library(grid)
library(gridExtra)
# =================== Section 0: Initialising ===================
setwd("C:/cs350/models/sprint1")
nsam <- 197
nrep <- 100
# get the 'wrong' data set, reassign, query model, predict on the spot
nattributes<-10
data_nuisance<-'rho80' # rho80 reflects average rho for drosophila
model_nuisance<-'rho40'
nuisance<-data_nuisance
load(paste('../../data/rdas/traindata_',nsam,'_',nrep,"_",nuisance,'.rda',sep=''))
nuisance<-model_nuisance
load(paste('../../data/rdas/model_boosting_nsam',nsam,"_",nuisance,'_nat',nattributes,'.rda',sep=''))
# quick fix
colnames(test)[2]<-'pop.1' # note here 
test<-test[,columns.model]
X_test<-test[,2:dim(test)[2]]
y_test<-test.Y
# load(paste('../../data/rdas/traindata_',nsam,'_',nrep,"_",nuisance,'.rda',sep=''))
# load(paste('../../data/rdas/prediction_analytic_nsam',nsam,"_",nuisance,'.rda',sep=''))
source('../functions_plot.R')

#pred_glm<-predict(m_glm,newdata=X_test)
#pred_glm<-summarisePred(pred_glm,y_test)
#p_glm<-plotPred(pred_glm,paste('glm:',model_nuisance,' model ',data_nuisance,' data',sep=''))
#p_glm

pred_gam<-predict(m_gam,newdata=X_test)
pred_gam<-summarisePred(pred_gam/5000,y_test/5000)
p_gam<-plotPred(pred_gam,paste('gam:',model_nuisance,' model ',data_nuisance,' data',sep=''))
p_gam

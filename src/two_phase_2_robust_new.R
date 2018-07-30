library(grid)
library(gridExtra)
library(caret)
library(mboost)
library(plyr)
library(ggplot2)
library(reshape2)
source('../functions_plot.R')
source('../functions_boosting.R')
setwd("C:/cs350/models/sprint1")

nsam <- 197
nrep <- 100

# first model and first test data
nuisance<-'pure'
load(paste('../../data/rdas/traindata_',nsam,'_',nrep,'_',nuisance,'.rda',sep=''))
# train the gam 
m_gam1<-traingam(train)
m_gam1_mus<-mus_train
m_gam1_sigmas<-sigmas_train
# obtain the corrected test data 
pure_test.X<-scaleback(test.X,mus_train,sigmas_train)
pure_test.X<-data.frame(pure_test.X)
pure_test.Y<-test.Y

# second model rho 40 
nuisance<-'rho40'
load(paste('../../data/rdas/traindata_',nsam,'_',nrep,'_',nuisance,'.rda',sep=''))
# train the gam 
m_gam2<-traingam(train)
m_gam2_mus<-mus_train
m_gam2_sigmas<-sigmas_train
# obtain the corrected test data 
rho40_test.X<-scaleback(test.X,mus_train,sigmas_train)
rho40_test.X<-data.frame(rho40_test.X)
rho40_test.Y<-test.Y

# rho 80
nuisance<-'rho80'
load(paste('../../data/rdas/traindata_',nsam,'_',nrep,'_',nuisance,'.rda',sep=''))
# train the gam 
m_gam3<-traingam(train)
m_gam3_mus<-mus_train
m_gam3_sigmas<-sigmas_train
# obtain the corrected test data 
rho80_test.X<-scaleback(test.X,mus_train,sigmas_train)
rho80_test.X<-data.frame(rho80_test.X)
rho80_test.Y<-test.Y

# rho 120
nuisance<-'rho120'
load(paste('../../data/rdas/traindata_',nsam,'_',nrep,'_',nuisance,'.rda',sep=''))
# train the gam 
m_gam4<-traingam(train)
m_gam4_mus<-mus_train
m_gam4_sigmas<-sigmas_train
# obtain the corrected test data 
rho120_test.X<-scaleback(test.X,mus_train,sigmas_train)
rho120_test.X<-data.frame(rho120_test.X)
rho120_test.Y<-test.Y

######### prediction ######### 
# pure gam ==> first row
# first scale using the model's scale
gam1_pure_test.X<-data.frame(scale(pure_test.X,center=m_gam1_mus,scale=m_gam1_sigmas))
gam1_rho40_test.X<-data.frame(scale(rho40_test.X,center=m_gam1_mus,scale=m_gam1_sigmas))
gam1_rho80_test.X<-data.frame(scale(rho80_test.X,center=m_gam1_mus,scale=m_gam1_sigmas))
gam1_rho120_test.X<-data.frame(scale(rho120_test.X,center=m_gam1_mus,scale=m_gam1_sigmas))

pred_gam1_pure<-predict(m_gam1,newdata=gam1_pure_test.X)
pred_gam1_pure<-summarisePred(pred_gam1_pure,pure_test.Y)

pred_gam1_rho40<-predict(m_gam1,newdata=gam1_rho40_test.X)
pred_gam1_rho40<-summarisePred(pred_gam1_rho40,rho40_test.Y)

pred_gam1_rho80<-predict(m_gam1,newdata=gam1_rho80_test.X)
pred_gam1_rho80<-summarisePred(pred_gam1_rho80,rho80_test.Y)

pred_gam1_rho120<-predict(m_gam1,newdata=gam1_rho120_test.X)
pred_gam1_rho120<-summarisePred(pred_gam1_rho120,rho120_test.Y)

# p1<-plotPred(pred_gam1_pure,'pure,pure',20,15)
p2<-plotPred2(pred_gam1_rho40,expression(paste('(a) ',rho,'=0 Boosted GAM, ', rho, '=40 test data',sep='')),20,15,0,225)
p3<-plotPred2(pred_gam1_rho80,expression(paste('(b) ',rho,'=0 Boosted GAM, ', rho, '=80 test data',sep='')),20,15,0,225)
p4<-plotPred2(pred_gam1_rho120,expression(paste('(c) ',rho,'=0 Boosted GAM, ', rho, '=120 test data',sep='')),20,15,0,225)

plots1<-grid.arrange(p2,p3,p4,nrow=1)

# rho40 gam ==> second row
# first scale using the model's scale
gam2_pure_test.X<-data.frame(scale(pure_test.X,center=m_gam2_mus,scale=m_gam2_sigmas))
gam2_rho40_test.X<-data.frame(scale(rho40_test.X,center=m_gam2_mus,scale=m_gam2_sigmas))
gam2_rho80_test.X<-data.frame(scale(rho80_test.X,center=m_gam2_mus,scale=m_gam2_sigmas))
gam2_rho120_test.X<-data.frame(scale(rho120_test.X,center=m_gam2_mus,scale=m_gam2_sigmas))

pred_gam2_pure<-predict(m_gam2,newdata=gam2_pure_test.X)
pred_gam2_pure<-summarisePred(pred_gam2_pure,pure_test.Y)

pred_gam2_rho40<-predict(m_gam2,newdata=gam2_rho40_test.X)
pred_gam2_rho40<-summarisePred(pred_gam2_rho40,rho40_test.Y)

pred_gam2_rho80<-predict(m_gam2,newdata=gam2_rho80_test.X)
pred_gam2_rho80<-summarisePred(pred_gam2_rho80,rho80_test.Y)

pred_gam2_rho120<-predict(m_gam2,newdata=gam2_rho120_test.X)
pred_gam2_rho120<-summarisePred(pred_gam2_rho120,rho120_test.Y)

# rho80 gam ==> third row
# first scale using the model's scale
gam3_pure_test.X<-data.frame(scale(pure_test.X,center=m_gam3_mus,scale=m_gam3_sigmas))
gam3_rho40_test.X<-data.frame(scale(rho40_test.X,center=m_gam3_mus,scale=m_gam3_sigmas))
gam3_rho80_test.X<-data.frame(scale(rho80_test.X,center=m_gam3_mus,scale=m_gam3_sigmas))
gam3_rho120_test.X<-data.frame(scale(rho120_test.X,center=m_gam3_mus,scale=m_gam3_sigmas))

pred_gam3_pure<-predict(m_gam3,newdata=gam3_pure_test.X)
pred_gam3_pure<-summarisePred(pred_gam3_pure,pure_test.Y)

pred_gam3_rho40<-predict(m_gam3,newdata=gam3_rho40_test.X)
pred_gam3_rho40<-summarisePred(pred_gam3_rho40,rho40_test.Y)

pred_gam3_rho80<-predict(m_gam3,newdata=gam3_rho80_test.X)
pred_gam3_rho80<-summarisePred(pred_gam3_rho80,rho80_test.Y)

pred_gam3_rho120<-predict(m_gam3,newdata=gam3_rho120_test.X)
pred_gam3_rho120<-summarisePred(pred_gam3_rho120,rho120_test.Y)



# rho120 gam ==> fourth row
# first scale using the model's scale
gam4_pure_test.X<-data.frame(scale(pure_test.X,center=m_gam4_mus,scale=m_gam4_sigmas))
gam4_rho40_test.X<-data.frame(scale(rho40_test.X,center=m_gam4_mus,scale=m_gam4_sigmas))
gam4_rho80_test.X<-data.frame(scale(rho80_test.X,center=m_gam4_mus,scale=m_gam4_sigmas))
gam4_rho120_test.X<-data.frame(scale(rho120_test.X,center=m_gam4_mus,scale=m_gam4_sigmas))

pred_gam4_pure<-predict(m_gam4,newdata=gam4_pure_test.X)
pred_gam4_pure<-summarisePred(pred_gam4_pure,pure_test.Y)

pred_gam4_rho40<-predict(m_gam4,newdata=gam4_rho40_test.X)
pred_gam4_rho40<-summarisePred(pred_gam4_rho40,rho40_test.Y)

pred_gam4_rho80<-predict(m_gam4,newdata=gam4_rho80_test.X)
pred_gam4_rho80<-summarisePred(pred_gam4_rho80,rho80_test.Y)

pred_gam4_rho120<-predict(m_gam4,newdata=gam4_rho120_test.X)
pred_gam4_rho120<-summarisePred(pred_gam4_rho120,rho120_test.Y)


# trained rho on rho test data 
u<-220
l<-0

p6<-plotPred2(pred_gam2_rho40,expression(paste('(a) ',rho,'=40 Boosted GAM, ', rho, '=40 test data',sep='')),20,15,l,u)
p7<-plotPred2(pred_gam2_rho80,expression(paste('(b) ',rho,'=40 Boosted GAM, ', rho, '=80 test data',sep='')),20,15,l,u)
p8<-plotPred2(pred_gam2_rho120,expression(paste('(c) ',rho,'=40 Boosted GAM, ', rho, '=120 test data',sep='')),20,15,l,u)

p10<-plotPred2(pred_gam3_rho40,expression(paste('(d) ',rho,'=80 Boosted GAM, ', rho, '=40 test data',sep='')),20,15,l,u)
p11<-plotPred2(pred_gam3_rho80,expression(paste('(e) ',rho,'=80 Boosted GAM, ', rho, '=80 test data',sep='')),20,15,l,u)
p12<-plotPred2(pred_gam3_rho120,expression(paste('(f) ',rho,'=80 Boosted GAM, ', rho, '=120 test data',sep='')),20,15,l,u)

p14<-plotPred2(pred_gam4_rho40,expression(paste('(g) ',rho,'=120 Boosted GAM, ', rho, '=40 test data',sep='')),20,15,l,u)
p15<-plotPred2(pred_gam4_rho80,expression(paste('(h) ',rho,'=120 Boosted GAM, ', rho, '=80 test data',sep='')),20,15,l,u)
p16<-plotPred2(pred_gam4_rho120,expression(paste('(i) ',rho,'=120 Boosted GAM, ', rho, '=120 test data',sep='')),20,15,l,u)

plots2<-grid.arrange(p6,p7,p8,p10,p11,p12,p14,p15,p16,nrow=3)

# trained rho on pure test data 
p5<-plotPred2(pred_gam2_pure,expression(paste('(a) ',rho,'=40 Boosted GAM, ', rho, '=0 test data',sep='')),20,15,0,270)
p9<-plotPred2(pred_gam3_pure,expression(paste('(b) ',rho,'=80 Boosted GAM, ', rho, '=0 test data',sep='')),20,15,0,270)
p13<-plotPred2(pred_gam4_pure,expression(paste('(c) ',rho,'=120 Boosted GAM, ', rho, '=0 test data',sep='')),20,15,0,270)

plots3 <- grid.arrange(p5,p9,p13,nrow=1)

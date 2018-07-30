library(plyr)
library(ggplot2)
library(grid)
library(gridExtra)
# =================== Section 0: Initialising ===================
setwd("C:/cs350/models/sprint1")
# change the input rda files to train different models and make different predictions
# this allows investigation into making prediction of datasets arising from different relationships of n and theta
nsam <- 100
nrep <- 100
nuisance<-'pure'
load(paste('../../data/rdas/traindata_',nsam,'_',nrep,'_',nuisance,'.rda',sep=''))
load('../../data/rdas/prediction_analytic_nsam100_pure_nat0.rda')
load('../../data/rdas/prediction_genetree.rda')
load('../../data/rdas/MOSTRECENTBOOST_100_100_pure.rda')
source('../functions_plot.R')

# ==================== Section 1: Evaluate Predictions of Boosting Models ====================
## gam
pred_gam <- predict(m_gam,newdata=test.X)
pred_gam <- summarisePred(pred_gam,test.Y)
limit <- max(pred_gam$mu+pred_gam$sd)
p1<-plotPred(pred_gam,"Boosted GAM",20,15)


# ==================== Section 2: Evaluate Predictions of Analytic Methods ====================
## Watterson
pred_w <- summarisePred(pred_w,test.Y)
p2<-plotPred2(pred_w,"(a) Watterson's Estimator",20,15,0,250)

## Tajima
pred_d <- summarisePred(pred_d,test.Y)
p3<-plotPred2(pred_d,"(b) Tajima's Estimator",20,15,0,250)

pred_fu <- summarisePred(pred_fu,test.Y)
p4<-plotPred2(pred_fu,"(c) Fu and Li's Estimator",20,15,0,250)


# ==================== Section 3: Evaluate Predictions of GeneTree ====================
load('../../data/rdas/prediction_genetree.rda') # to edit
# change the following to melt then subject it to the same summarisePred as the rest
pred_gt<-melt(data,id='seed')
pred_gt$true<-as.numeric(as.character(substr(pred_gt$variable,2,length(pred_gt$variable))))
pred_gt<-pred_gt[,c(3,4)]
colnames(pred_gt)<-c('hat','true')
# routine
pred_gt<-summarisePred(pred_gt[,1],pred_gt[,2])
p5<-plotPred2(pred_gt,"(d) GENETREE",20,15,0,250)


plots<-grid.arrange(p2,p3,p4,p5,nrow=2,ncol=2)

# ==================== Section 4: RMSE Comparison ====================
library(reshape2)
mse <- data.frame(list(true=pred_gam$true,
                        GAM=pred_gam$mse,
                        W=pred_w$mse,
                        T=pred_d$mse,
                        F=pred_fu$mse,
                        GT=pred_gt$mse
))
mse<-melt(mse,id='true')
colnames(mse)[2]<-"Methods"
c1<-ggplot(mse)+
  geom_line(aes(x=true,y=value,colour=Methods))+
  scale_colour_manual(values=c("red","green","blue","black","purple","brown","orange","pink","yellow"))+
  ggtitle("Comparison of Point-wise MSE")+
  theme(plot.title=element_text(size=15,hjust = 0.5)) +
  labs(x=expression(paste("True ",theta,sep='')),
       y="MSE")

# ==================== Section 5: Bias Comparison ====================
## Bias 
bias <- data.frame(list(true=pred_gam$true,
                        GAM=pred_gam$bias,
                        W=pred_w$bias,
                        D=pred_d$bias,
                        F=pred_fu$bias,
                        GT=pred_gt$bias
))
bias<-melt(bias,id='true')
colnames(bias)[2]<-"Models"
c2<-ggplot(bias)+
  geom_line(aes(x=true,y=value,colour=Models))+
  scale_colour_manual(values=c("red","green","blue","black","purple","brown","orange","pink"))+
  ggtitle("Comparison of Point-wise Bias")+
  theme(plot.title=element_text(size=15,hjust = 0.5)) +
  labs(x=expression(paste("True ",theta,sep='')),
       y="Bias")

## Bias^2
bias2 <- data.frame(list(true=pred_gam$true,
                        GAM=pred_gam$bias^2,
                        W=pred_w$bias^2,
                        T=pred_d$bias^2,
                        F=pred_fu$bias^2,
                        GT=pred_gt$bias^2
))
bias2<-melt(bias2,id='true')
colnames(bias2)[2]<-"Methods"
c3<-ggplot(bias2)+
  geom_line(aes(x=true,y=value,colour=Methods))+
  scale_colour_manual(values=c("red","green","blue","black","purple","brown","orange"))+
  ggtitle("Comparison of Point-wise Bias Squared")+
  theme(plot.title=element_text(size=15,hjust = 0.5)) +
  labs(x=expression(paste("True ",theta,sep='')),
       y="Bias squared")

## variance
var <- data.frame(list(true=pred_gam$true,
                        GAM=pred_gam$var,
                        W=pred_w$var,
                        T=pred_d$var,
                        F=pred_fu$var,
                        GT=pred_gt$var
))
var<-melt(var,id='true')
colnames(var)[2]<-"Methods"
c4<-ggplot(var)+
  geom_line(aes(x=true,y=value,colour=Methods))+
  scale_colour_manual(values=c("red","green","blue","black","purple","brown","orange"))+
  ggtitle("Comparison of Point-wise Variance")+
  theme(plot.title=element_text(size=15,hjust = 0.5)) +
  labs(x=expression(paste("True ",theta,sep='')),
       y="Variance")


# ==================== Saving Sections ====================
save(p1,p2,p3,p4,p5,p6,p7,p8,c1,c2,c3,c4,c5,file=paste('../../data/rdas/prediction_evaluation_graphs_nsam',nsam,'.rda',sep=''))
save(pred_glm,pred_gam1,pred_gam2,pred_w,pred_d,pred_fu,pred_gt,file=paste('../../data/rdas/prediction_summarised_nsam',nsam,'.rda',sep=''))


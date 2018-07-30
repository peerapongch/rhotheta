library(plyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)
source('../functions_plot.R')
nsam <- 100
nrep <- 100
nuisance<-'pure'
load(paste('../../data/rdas/traindata_',nsam,'_',nrep,'_',nuisance,'.rda',sep=''))

# tree
m_tre<-blackboost(theta~.,
                  data=train,
                  family=Gaussian(),
                  control=boost_control(mstop=500,
                                        nu=0.1,
                                        trace=TRUE
                  )
)
system.time(
  cvm <- cvrisk(m_tre, folds = cv(weights = model.weights(m_tre),
                                  type = "kfold"))
)
m_tre[mstop(cvm)]

# gam
m_gam <- gamboost(theta~.,
                  data=train,
                  control=boost_control(mstop=500,
                                        trace=TRUE
                  ),
                  baselearner = "bbs",
                  family=Gaussian(),
                  dfbase=4
)
system.time(
  cvm1 <- cvrisk(m_gam, folds = cv(weights = model.weights(m_gam),
                                   type = "kfold"))
)
m_gam[mstop(cvm1)]

# glm
# m1: glm fitting
m_glm <- glmboost(theta~.,
                  data=train,
                  family=Gaussian(), #what about GaussReg()?
                  control = boost_control(mstop = 500,
                                          nu = 0.1
                  )
)

system.time(
  cvm1 <- cvrisk(m_glm, folds = cv(weights = model.weights(m_glm),
                                   type = "kfold"))
)

m_glm[mstop(cvm1)]


# save(m_glm,m_gam,m_tre,file='../../data/rdas/MOSTRECENTBOOST_100_100_pure.rda')
# load('../../data/rdas/MOSTRECENTBOOST_100_100_pure.rda')
# now that we have the models 

# output the lins graphs
pred_tre<-predict(m_tre,newdata=test.X)
pred_gam<-predict(m_gam,newdata=test.X)
pred_glm<-predict(m_glm,newdata=test.X)

pred_tre<-summarisePred(pred_tre,test.Y)
pred_gam<-summarisePred(pred_gam,test.Y)
pred_glm<-summarisePred(pred_glm,test.Y)


p1<-plotPred(pred_glm,'(a) Boosted GLM',20,15)
p2<-plotPred(pred_gam,'(b) Boosted GAM',20,15)
p3<-plotPred(pred_tre,'(c) Boosted TRE',20,15)

plots<-grid.arrange(p1,p2,p3,nrow=1,ncol=3)

# use ggplot2 to make the fancy pointwise graph 

# mse 
mse <- data.frame(list(true=pred_glm$true,
                        GLM=pred_glm$mse,
                        GAM=pred_gam$mse,
                        TRE=pred_tre$mse
))
mse<-melt(mse,id='true')
colnames(mse)[2]<-"Models"
c1<-ggplot(mse)+
  geom_line(aes(x=true,y=value,colour=Models))+
  scale_colour_manual(values=c("red","blue","green"))+
  ggtitle("Point-wise MSE")+
  theme(plot.title=element_text(size=15,hjust = 0.5)) +
  labs(x=expression(paste("True ",theta,sep='')),
       y="MSE") +
  geom_vline(xintercept=100, linetype="dotted")
  

# bias
bias <- data.frame(list(true=pred_glm$true,
                        GLM=pred_glm$bias,
                        GAM=pred_gam$bias,
                        TRE=pred_tre$bias
))
bias<-melt(bias,id='true')
colnames(bias)[2]<-"Models"
c2<-ggplot(bias)+
  geom_line(aes(x=true,y=value,colour=Models))+
  scale_colour_manual(values=c("red","blue","green"))+
  ggtitle("Point-wise Bias")+
  theme(plot.title=element_text(size=15,hjust = 0.5)) +
  labs(x=expression(paste("True ",theta,sep='')),
       y="Bias") +
  geom_vline(xintercept=100, linetype="dotted")

# var
var <- data.frame(list(true=pred_glm$true,
                        GLM=pred_glm$var,
                        GAM=pred_gam$var,
                        TRE=pred_tre$var
))
var<-melt(var,id='true')
colnames(var)[2]<-"Models"
c3<-ggplot(var)+
  geom_line(aes(x=true,y=value,colour=Models))+
  scale_colour_manual(values=c("red","blue","green"))+
  ggtitle("Point-wise Variance")+
  theme(plot.title=element_text(size=15,hjust = 0.5)) +
  labs(x=expression(paste("True ",theta,sep='')),
       y="Variance") +
  geom_vline(xintercept=100, linetype="dotted")

bias2 <- data.frame(list(true=pred_glm$true,
                        GLM=pred_glm$bias^2,
                        GAM=pred_gam$bias^2,
                        TRE=pred_tre$bias^2
))
bias2<-melt(bias2,id='true')
colnames(bias2)[2]<-"Models"
c4<-ggplot(bias2)+
  geom_line(aes(x=true,y=value,colour=Models))+
  scale_colour_manual(values=c("red","blue","green"))+
  ggtitle("Point-wise Bias Squared")+
  theme(plot.title=element_text(size=15,hjust = 0.5)) +
  labs(x=expression(paste("True ",theta,sep='')),
       y="Bias Squared") +
  geom_vline(xintercept=100, linetype="dotted")

plots2 <- grid.arrange(c1,c2,c4,c3,nrow=2,ncol=2)

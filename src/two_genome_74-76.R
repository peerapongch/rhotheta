library(mboost)
library(reshape2)
library(ggplot2)
source('../functions_boosting.R')
source('../functions_plot.R')
pop<-'ZI'
tid<-'2L'
width<-'5000'

# first 7.4-7.5
frompos<-'7400000'
topos<-'7500000'
nuisance<-'rho387'
filename<-paste('../../data/genomes/',pop,'_Chr',tid,'_',frompos,'_',topos,'_window',width,'.csv',sep='')
first_data<-read.csv(filename)
# normalise by the training mean and sd
load(paste('../../data/rdas/traindata_197_100_',nuisance,'.rda',sep=''))
first_test.X<-first_data[,columns.X]
first_test.X<-data.frame(scale(first_test.X,center=mus_train,scale=sigmas_train))
# train the boosting model for this locus 
first_gam<-traingam(train)
pred_first_gam<-predict(first_gam,newdata=first_test.X)

# then 7.5-7.6
frompos<-'7500000'
topos<-'7600000'
nuisance<-'rho105'
filename<-paste('../../data/genomes/',pop,'_Chr',tid,'_',frompos,'_',topos,'_window',width,'.csv',sep='')
second_data<-read.csv(filename)
# standardise
load(paste('../../data/rdas/traindata_197_100_',nuisance,'.rda',sep=''))
second_test.X<-second_data[,columns.X]
second_test.X<-data.frame(scale(second_test.X,center=mus_train,scale=sigmas_train))
# boosting model
second_gam<-traingam(train)
pred_second_gam<-predict(second_gam,newdata=second_test.X)

##### merging the two parts and plot ######
positions<-cumsum(c(first_data$bp,second_data$bp))+7400000
pred_gam<-c(pred_first_gam,pred_second_gam)/5000
pred_w<-c(first_data$theta_watterson,second_data$theta_watterson)/5000
pred_t<-c(first_data$theta_tajima,second_data$theta_tajima)/5000
pred_f<-c(first_data$theta_fu,second_data$theta_fu)/5000

pred <- data.frame(list(region=positions,
                        #glm=pred_glm,
                        GAM=pred_gam,
                        W=pred_w,
                        T=pred_t,
                        F=pred_f
                        )
                   )
pred<-melt(pred,id='region')
colnames(pred)[2]<-"Methods"
locustheta<-ggplot(pred)+
  geom_line(aes(x=region,y=value,colour=Methods))+
  scale_colour_manual(values=c("red","blue","purple","black"))+
  theme(plot.title=element_text(size=20,hjust = 0.5),
        text=element_text(size=15)) +
  ggtitle(expression(paste('Estimates of ',theta,' per base pair along 7.4-7.6 Mb stretch at 5 kb intervals',sep=''))
  ) +
  labs(
    x="bp",
    y=expression(paste("Estimated ",theta,sep=''))
  )+
  geom_hline(yintercept=0.016, linetype="dashed")
locustheta

plot(seq(7405000,7600000,5000),c(first_data$test_tajima_d,second_data$test_tajima_d), type='l')
abline(h=-2)

#### plotting separately ####
## 7.4-7.5
llim<-7400000
ulim<-7500000
positions <- cumsum(first_data$bp)+llim-2500
pred <- data.frame(list(region=positions,
                        #glm=pred_glm,
                        GAM=pred_first_gam/5000,
                        W=first_data$theta_watterson/5000,
                        T=first_data$theta_tajima/5000,
                        F=first_data$theta_fu/5000
)
)
pred<-melt(pred,id='region')
colnames(pred)[2]<-"Methods"
locustheta1<-ggplot(pred)+
  geom_line(aes(x=region,y=value,colour=Methods))+
  scale_colour_manual(values=c("red","blue","purple","black"))+
  theme(plot.title=element_text(size=20,hjust = 0.5),
        text=element_text(size=15)) +
  ggtitle(expression(paste('Estimates of ',theta,' per base pair along 7.4-7.5 Mb stretch at 5 kb intervals',sep=''))
  ) +
  labs(
    x="bp",
    y=expression(paste("Estimated ",theta,sep=''))
  )+
  geom_hline(yintercept=0.016, linetype="dashed")+
  xlim(llim,ulim)
locustheta1

## 7.5-7.6
llim<-7500000
ulim<-7600000
positions <- cumsum(first_data$bp)+llim-2500
pred <- data.frame(list(region=positions,
                        #glm=pred_glm,
                        GAM=pred_second_gam/5000,
                        W=second_data$theta_watterson/5000,
                        T=second_data$theta_tajima/5000,
                        F=second_data$theta_fu/5000
)
)
pred<-melt(pred,id='region')
colnames(pred)[2]<-"Methods"
locustheta2<-ggplot(pred)+
  geom_line(aes(x=region,y=value,colour=Methods))+
  scale_colour_manual(values=c("red","blue","purple","black"))+
  theme(plot.title=element_text(size=20,hjust = 0.5),
        text=element_text(size=15)) +
  ggtitle(expression(paste('Estimates of ',theta,' per base pair along 7.5-7.6 Mb stretch at 5 kb intervals',sep=''))
  ) +
  labs(
    x="bp",
    y=expression(paste("Estimated ",theta,sep=''))
  )+
  geom_hline(yintercept=0.016, linetype="dashed")+
  xlim(llim,ulim)
locustheta2

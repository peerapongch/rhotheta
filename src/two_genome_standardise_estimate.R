# Note: crude selection of first 50 etas used to test the extensibility of pretrained models on arbitrary nsam 
# ideally, new data with nsam=n_genome should be generated and new models should be trained
library(mboost)
library(reshape2)
library(ggplot2)
source('../functions_boosting.R')
source('../functions_plot.R')

# === Initialise ====
pop<-'ZI'
tid<-'2L'
frompos<-'7400000'
topos<-'7600000'
width<-'5000'
filename<-paste('../../data/genomes/',pop,'_Chr',tid,'_',frompos,'_',topos,'_window',width,'.csv',sep='')
nuisance<-'pure'

# === New Data ====
new_data<-read.csv(filename)
# normalise by the training mean and sd
load(paste('../../data/rdas/traindata_197_100_',nuisance,'.rda',sep=''))
# subset the covariates 
X<-new_data[,columns.X]
X<-t(apply(X,FUN=function(x){(x-mus_train)/sigmas_train},MARGIN=1))
X<-as.data.frame(X)

# === estimation per bp ====
load(paste('../../data/rdas/model_boosting_nsam197_',nuisance,'_nat0.rda',sep=''))
# includes model specific subset of columns (generalised from testing for robustness)
X<-X[,columns.model[2:length(columns.model)]]
# with analytic 
pred_w<-new_data$theta_watterson/new_data$bp
pred_d<-new_data$theta_tajima/new_data$bp
pred_fu<-new_data$theta_fu/new_data$bp # we can choose freely either fu or eta1 here 

# prediction with glm
#pred_glm<-predict(m_glm,newdata=X)/new_data$bp

# prediction with gam 
pred_gam<-predict(m_gam,newdata=X)/new_data$bp


# === visualisation 1: boosting prediction against analytic ====
# just a quick one 
#plot(pred_fu,pred_gam,xlim=c(0,0.05),ylim=c(0,0.05))
#abline(a=0,b=1)

# === visualisation 2: plot of theta per bp along chromosome ====
# source('functions_plot.R')
pred <- data.frame(list(region=as.numeric(row.names(new_data))*as.numeric(width)+as.numeric(frompos),
                        #glm=pred_glm,
                        GAM=pred_gam,
                        W=pred_w,
                        T=pred_d,
                        F=pred_fu
                        )
                   )
pred<-melt(pred,id='region')
colnames(pred)[2]<-"Methods"
locustheta<-ggplot(pred)+
  geom_line(aes(x=region,y=value,colour=Methods))+
  scale_colour_manual(values=c("red","blue","purple","black"))+
  theme(plot.title=element_text(size=20,hjust = 0.5),
        text=element_text(size=15)) +
  ggtitle(expression(paste('Estimates of ',theta,' per base pair along 7.4-7.6 Mb at 5 kb intervals',sep=''))
  ) +
  labs(
    x="bp",
    y=expression(paste("Estimated ",theta,sep=''))
    )+
  geom_hline(yintercept=0.016, linetype="dashed")
locustheta


#plot(new_data$test_tajima_d, type='l')
#plot(new_data$test_fu_li_d, type='l')

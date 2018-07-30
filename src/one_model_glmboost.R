mse <- function(ytrue,yhat){
  return(mean((ytrue-yhat)^2))
}

# implements the boosted glm from the summarised data

load('../../data/repo2/data_100_100.rda')

# set the plot pane
# par(mfrow=c(1,1))

# m1: glm fitting
m_glm <- glmboost(theta~.,
               data=train,
               family=Gaussian(), #what about GaussReg()?
               control = boost_control(mstop = 500,
                                       nu = 0.1
                                       )
               )
pred_glm <- predict(m_glm, newdata=test.X)
plot(test.Y, theta_hat_m1, xlim=c(-50,200), ylim=c(-50,200),asp=1)
abline(a=0,b=1)
# explore the risk evaluation

# by CV
system.time(
  cvm1 <- cvrisk(m_glm, folds = cv(weights = model.weights(m_glm),
                                 type = "kfold"))
)

m_glm[mstop(cvm1)]

mse(predict(m_glm[mstop(cvm1)],newdata=test.X),test.Y)

# by AIC
system.time(
  aic_glm<-AIC(m_glm)
)

mstop(aic_glm)

mse(predict(m_glm[mstop(aic_glm)],newdata=test.X),test.Y)

cvm <- cvrisk(m1)
plot(cvm)
m1[mstop(cvm)]
plot(m1)

### coef paths of selected coefficients
m_glm[142]
library(ggplot2)
library(reshape2)
library(directlabels)
cf1 <- coef(m_glm, which = c(2,3,6,7,8,9), aggregate = "cumsum")
tmp <- sapply(cf1, function(x) x)
#matplot(tmp, type = "l", main = "Coefficient Paths")
colnames(tmp) <- c('S','K','eta1','eta2','eta3','eta4')
coef <- melt(tmp)
colnames(coef)[2]<-"Variable"
#coef<-coef[coef$value!=0,]
c1<-ggplot(coef)+
  geom_line(aes(x=Var1,y=value,colour=Variable))+
  scale_colour_manual(values=c("black","black","black","black","black","black"))+
  ggtitle("Coefficient Paths")+
  theme(plot.title=element_text(size=13,hjust = 0.5)) +
  labs(y=expression(paste("Coefficient estimate",sep='')),
       x="Number of boosting iteration")+
  xlim(0,150)+
  ylim(0,18)
direct.label(c1, list(last.points, hjust = -0.3, vjust = 0.5))




# dont care from here on

# m1_cv: glm with 10-fold cross validation 
fitControl <- trainControl(method="cv",
                           number=10
)
system.time(
  {
    m1_cv <- train(theta~.,
               data=train,
               method="glmboost",
               trControl=fitControl,
               tuneGrid=expand.grid(mstop = seq(10,500,10), prune = c("yes","no"))
               )
  }
)
trellis.par.set(caretTheme())
plot(m1_cv)
# why no difference between pruning and non-pruning methods?
m1_cv_opt <- m1_cv$finalModel
theta_hat_m1_cv <- predict(m1_cv_opt,newdata=test.X)
plot(test.Y, theta_hat_m1_cv,xlim=c(0,150), ylim=c(-100,500))
abline(a=0,b=1)

# evaluation
rmse_m0 <- sqrt(mean((theta_hat_m0-test.Y)^2))
rmse_m1 <- sqrt(mean((theta_hat_m1-test.Y)^2))
rmse_m1_cv <- sqrt(mean((theta_hat_m1_cv-test.Y)^2))

# coefficients
data.frame("glm"=coef(m1),"glm_cv"=coef(m1_cv_opt))

# evaluation 
n <- dim(train.X)[1]
res <- seq(400,450,1)
aics <- rep(-1,length(res))
aics_c <- rep(-1,length(res))
bics <- rep(-1,length(res))
ps <- rep(-1,length(res))
ll[i] <- rep(-1,length(res))
for(i in 1:length(res)){
  m1 <- glmboost(theta~.,
                 data=train,
                 family=Gaussian(), 
                 control = boost_control(mstop = res[i],
                                         nu = 0.1
                                         )
                 )
  aics[i] <- AIC(m1)[1]
  aics_c[i] <- AIC(m1, method="corrected")
  bics[i] <- log(n)*length(coef(m1))-2*logLik(m1)
  ps[i] <- length(coef(m1))
  rss[i] <- sqrt(mean(resid(m1)^2))
  ll[i] <- logLik(m1)
}
# something weird
aic_diy <- 2*(ps)-2*ll
plot(res,aic_diy)
plot(res,2*ps+n*log(rss/n)-2*c) 
plot(res,aics)

par(mar = c(5,5,2,5))
plot(res,aics,"l",xlab = 'm',ylab='AIC',ylim=c(6,6.5))
lines(x=res,y=aics_c,col='red')
par(new=T)
plot(x=res,y=bics,col='blue',axes=F, xlab=NA, ylab=NA, type='l')
axis(side = 4)
mtext("BIC", side = 4, line = 3)
save(res,aics,aics_c,bics,file="ics.rda")


# why does AIC flatline after iteration 409?
aic <- AIC(m1)[1]
ll <- (aic-2*length(coef(m1)))/2
rss <- sqrt(mean(resid(m1)^2))
n <- dim(train.X)[1]
c <- ll + n/2*log(rss/n)


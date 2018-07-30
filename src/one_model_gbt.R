# gradient boosted tree with mboost 
library(mboost)
# =================== Section 0: Initialising ===================
setwd("C:/cs350/models/sprint1")
# change the input rda files to train different models and make different predictions
# this allows investigation into making prediction of datasets arising from different relationships of n and theta
nsam <- 100
nrep <- 100
nuisance <- 'pure'
load(paste('../../data/rdas/traindata_',nsam,'_',nrep,'_','pure','.rda',sep=''))

# with default ctree settings
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

pred_tre<-predict(m_tre,newdata=test.X)
plot(test.Y,pred_tre,xlim=c(0,150),ylim=c(0,150))
abline(a=0,b=1)


mse(test.Y,pred_tre)
sqrt(mean((pred_tre-test.Y)^2))

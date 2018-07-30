library(caret)
library(mboost)
library(plyr)
library(ggplot2)
# =================== Section 0: Initialising ===================
setwd("C:/cs350/models/sprint1")
# change the input rda files to train different models and make different predictions
# this allows investigation into making prediction of datasets arising from different relationships of n and theta
nsam <- 100
nrep <- 100
load(paste('../../data/repo2/data_',nsam,'_',nrep,'.rda',sep=''))

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

plot(cvm1)

mse(predict(m_gam[mstop(cvm1)],newdata=test.X),test.Y)

system.time(
aic<-AIC(m_gam)
)

mse(predict(m_gam[mstop(aic)],newdata=test.X),test.Y)


## shortcuts?
df<-scale(train.X,center=FALSE,scale=1/sigmas_train)
df<-scale(df,center=-mus_train,scale=FALSE)
temp_train<-data.frame(cbind(train[,1],df))
colnames(temp_train)[1]<-'theta'

m_gam <- gamboost(theta~.,
                  data=temp_train,
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

# produce splines plots
library(grid)
library(gridExtra)

par(cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
par(mfrow=c(3,2))
par(mar=c(5,7,4,2))
plot(m_gam,which=1,type='l',xlab = 'S')
plot(m_gam,which=2,type='l',xlab = 'K')
plot(m_gam,which=3,type='l',xlab = 'Pi')
plot(m_gam,which=5,type='l',xlab = 'eta1')
plot(m_gam,which=6,type='l',xlab = 'eta2')
plot(m_gam,which=7,type='l',xlab = 'eta3')









# dont care about below 

# =================== Section 2: some slight changes of knots ===================
# say we only focus on two covariates S_n and H_n 
# m2 benchmark 
m2 <- gamboost(theta~.,
               data=train[,1:4],
               control=boost_control(mstop=500,
                                     trace=TRUE
               ),
               baselearner = "bbs",
               family=Gaussian(),
               dfbase=4
)
cvm2 <- cvrisk(m2, folds = cv(weights = model.weights(m2),
                                 type = "kfold"))
plot(cvm2)
m2[mstop(cvm2)]
plot(m2)

# m3 with prespecified splines 
x<-train[,2]
s<-bbs(x,
       boundary.knots=c(min(x),max(x)),
       knots=length(x)
       )
h<-bbs(train[,3],
       boundary.knots=c(min(train[,3]),max(train[,3])),
       knots=length(train[,3])
)
p<-bbs(train[,4],
       boundary.knots=c(min(train[,4]),max(train[,4])),
       knots=length(train[,4])
)
m3<-gamboost(theta~s+h+p,
             data=train,
             control=boost_control(mstop=500,
                                   trace=TRUE,
                                   nu=0.01
                                   ),
             family=Gaussian()
             )
cvm3 <- cvrisk(m3, folds = cv(weights = model.weights(m3),
                              type = "kfold"))
plot(cvm3)
m2[mstop(cvm3)]
plot(m3)

m4<-gamboost(theta~s,
             data=train,
             control=boost_control(mstop=500,
                                   trace=TRUE,
                                   nu=0.001
             ),
             family=Gaussian()
)
plot(m4)
# prediction
pred_m2 <- predict(m2, newdata = test.X)
pred_m4 <- predict(m4, newdata = list(x=test[,2]))
# the beyond boundary knots stuff 
x2<-test.X[,2]
sum(x2[x2>max(x)])
sum(x2<min(x))

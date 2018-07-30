# implement alternative glm to investigate the effect of varying the ordering of covariates in the boosting algorithm
load('../../data/repo2/data_100_100.rda')
library(mboost)

# 1. reverse the ordering of columns 
new_train.X<-train.X[,seq(dim(train.X)[2],1)]
new_train<-cbind(train.Y,new_train.X)
colnames(new_train)[1]<-'theta'
# fitting the glmboost with 10-fold cross validation
m1 <- glmboost(theta~.,
               data=new_train,
               family=Gaussian(), #what about GaussReg()?
               control = boost_control(mstop = 500,
                                       nu = 0.1
               )
)
summary(m1)
plot(m1)
# the result confirms that the boosted model is invariant to changes in the ordering of covariates (i.e. the ordering of variables entering the model)

# 2. removing covaritate highly correlated with H_diversity
cov(train[,2:11])
# identified as H_n
index<-which(colnames(new_train)=='H_n')
new_train2<-new_train[,-index]
#boosting glm
m2 <- glmboost(theta~.,
               data=new_train2,
               family=Gaussian(), #what about GaussReg()?
               control = boost_control(mstop = 500,
                                       nu = 0.1
               )
)
summary(m2)
plot(m2)
# result in greater significance of S_n and X1 and sligth increase in significance of H_diversity, the magnitude is still much lower than H_n 
# therefore, can pretty much disregard H_diversity
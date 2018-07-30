mse <- function(ytrue,yhat){
  return(mean((ytrue-yhat)^2))
}
scaleback <- function(data,mus,sigmas){
  df<-scale(data,center=FALSE,scale=1/sigmas)
  df<-scale(df,center=-mus,scale=FALSE)
  return(df)
}
traingam<- function(train){
  m_gam1 <- gamboost(theta~.,
                     data=train,
                     control=boost_control(mstop=500,
                                           trace=TRUE
                     ),
                     baselearner = "bbs",
                     family=Gaussian(),
                     dfbase=4
  )
  system.time(
    cvm1 <- cvrisk(m_gam1, folds = cv(weights = model.weights(m_gam1),
                                      type = "kfold"))
  )
  m_gam1[mstop(cvm1)]
  return(m_gam1)
}
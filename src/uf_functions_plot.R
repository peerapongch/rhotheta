library(plyr)
library(ggplot2)

summarisePred<-function(pred,true){
  # df must have true and hat columns
  df<-data.frame(list(hat=pred,true=true))
  df <- ddply(df,
              "true", # this the groupby column
              summarise,
              N=length(hat),
              mu=mean(hat),
              sd=sd(hat),
              # se=sd/sqrt(N),
              mse=mean((true-hat)^2),
              var=sd^2,
              bias=mean(hat-true),
              # bias2=bias^2,
              # zero=mse-bias2-var,
              rmse=sqrt(mse)
              )
  return(df)
}

plotPred<-function(data,title,title_size,element_size){
  u_limit <- max(data$mu+data$sd,data$true)
  l_limit <- min(data$mu-data$sd,data$true,c(0))
  return(
    ggplot(data, aes(x=true, y=mu)) + 
    # summary of prediction
    geom_errorbar(aes(ymin=mu-sd, ymax=mu+sd), width=.1) +
    # predicted against true
    geom_abline(intercept=0,slope=1) + 
    geom_point() +
    theme(plot.title=element_text(size=title_size,hjust = 0.5),
      text=element_text(size=element_size)) +
    ggtitle(title) +
    labs(x=expression(paste("True ",theta,sep='')),y=expression(paste("Estimated ",theta,sep=''))) +
    coord_fixed()+
    xlim(l_limit,u_limit)+
    ylim(l_limit,u_limit)
  )
}

plotPred2<-function(data,title,title_size,element_size,llim,ulim){
  u_limit <- ulim
  l_limit <- llim
  return(
    ggplot(data, aes(x=true, y=mu)) + 
    # summary of prediction
    geom_errorbar(aes(ymin=mu-sd, ymax=mu+sd), width=.1) +
    # predicted against true
    geom_abline(intercept=0,slope=1) + 
    geom_point() +
    theme(plot.title=element_text(size=title_size,hjust = 0.5),
      text=element_text(size=element_size)) +
    ggtitle(title) +
    labs(x=expression(paste("True ",theta,sep='')),y=expression(paste("Estimated ",theta,sep=''))) +
    coord_fixed()+
    xlim(l_limit,u_limit)+
    ylim(l_limit,u_limit)
  )
}

# plotLocus<-function(pred_w,pred_d,pred_fu,pred_glm,pred_gam,pred_tre,pred_gt){

# }
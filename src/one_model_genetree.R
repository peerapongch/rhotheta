# this script executes prediction with genetree from start to end 

# ==================== Set Parameters ===================
# for simulation 
nrepeat=100
nsam=200
thetas <- c(seq(10,150,10))
setwd("C:/cs350/models/sprint1")
pathtogtdata<-"../../data/training_pure/gt/"
# for genetree
library(plyr)
mutation=10
delta<-10 # this is to cheat a little bit and to help the code run faster and guide genetree to make a better estimation
surface_points<-20
pathtogenetree="../../genetree90/genetree "
run=100000
seed=0 # arbitrary value, resetted at every iteration
repeatset<-seq(1,nrepeat,1)
data<-data.frame(list(seed=repeatset))

# ==================== Do: ===================
sink(file="../gt_sink.log")
system.time(
  for(theta in thetas){
    mles<-c()
    for (rep in repeatset){
      #print(paste("ITERATION: ",rep,sep=""))
      seed<-rep
      # ==================== Simulate Data from ms ===================
      #print(paste("simulating for theta = ",theta,sep=""))
      filename<-paste(pathtogtdata,
                      "gt",
                      nsam,
                      "_t",
                      theta,
                      ".txt",
                      sep=""
      )
      cmd<-paste("bash -c \" ../../msdir/ms ",
                 nsam,
                 " ",
                 1,
                 " -seeds ",
                 seed,
                 " ",
                 seed,
                 " ",
                 seed,
                 " -t ",
                 theta,
                 " | tail -n +7 > ",
                 filename,
                 "\"",
                 sep=""
      )
      system(cmd)
      
      # ==================== Reformat and Run genetree ===================
      
      # reformatting the data 
      #print(paste("now theta is ",theta, sep=""))
      filename<-paste(pathtogtdata,"gt",nsam,"_t",theta,".txt",sep="")
      df<-read.csv(filename,header=FALSE,colClasses="character")
      df<-ddply(df,
                'V1',
                summarise,
                n=length(V1))
      df['colon']=rep(':',dim(df)[1])
      df<-df[,c(2,3,1)]
      write.table(df$V1,file="temp.csv",sep=" ",row.names=FALSE,col.names=FALSE)
      cmd <- "bash -c \"sed s/\\\\\\\"//g temp.csv > temp.dat\""
      system(cmd)
      cmd <- "bash -c \"sed -e 's/\\(.\\)/\\1 /g' temp.dat > temp.txt\""
      system(cmd)
      df2<-read.table("temp.txt",sep=" ",fill=TRUE)
      df2 <- df2[rowSums(is.na(df2))<ncol(df2),]
      df2 <- df2[,colSums(is.na(df2))<nrow(df2)]
      result<-cbind(df[,c(1,2)],df2)
      writetofile=paste(pathtogtdata,"gt",nsam,"_t",theta,".csv",sep="")
      write.table(result,file=writetofile,sep=" ",col.names=FALSE,row.names=FALSE)
      writetofilefinal=paste(pathtogtdata,"gt",nsam,"_t",theta,".dat",sep="")
      cmd <- paste("bash -c \"sed s/\\\\\\\"//g ", writetofile," > ", writetofilefinal, "\"",sep="")
      system(cmd)
      
      # running genetree
      maxtype<-dim(df2)[1]
      #print(paste("now genetreeing for theta ",theta,"; maxtype is ",maxtype,sep=""))
      # generate tree file and run genetree
      treefile<-paste(pathtogtdata,
                      "tree_",
                      nsam,
                      "_t",
                      theta,
                      ".dat",
                      sep="")
      surffile<-paste(pathtogtdata,
                      "surf_",
                      nsam,
                      "_t",
                      theta,
                      ".dat",
                      sep="")
      incidencefile=paste(pathtogtdata,"gt",nsam,"_t",theta,".dat",sep="")
      cmd<-paste("bash -c \"",
                 pathtogenetree,
                 treefile,
                 " ",
                 mutation,
                 " ",
                 run,
                 " ",
                 seed,
                 " -j ",
                 incidencefile,
                 " -x 3000",
                 " -y ",
                 maxtype,
                 " -f ",
                 surffile,
                 " -g ",
                 theta-delta+1,
                 " ",
                 theta+delta
                 ," ",
                 surface_points
                 #," -2 \""
                 ," \""
      )
      system(cmd)
      
      # open surffile and figure out the theta_mle
      surf<-read.table(surffile,sep=" ",header=FALSE)
      ml_index<-which.max(surf$V4)
      mle<-surf$V3[ml_index]
      print(paste(rep," : ",mle,sep=""))
      # store values
      mles<-c(mles,mle)
      print(mles)
    }
    # store in dataframe 
    name<-paste('t',theta,sep="")
    data[name]<-mles
  }
)
save(data,file=paste('../../data/rdas/prediction_genetree_'+nsam+'.rda',sep=''))
sink()
# ==================== End ===================

# this script takes in the input from Hudson's ms and comput summary statistics before saving the data into a csv file
# from simulate_ms.R
# to data.R
compute_ss <- function(nsam,nrep,theta_lower,theta_upper,increment,repo){
  # this is the simple serialised version
  # further modification may be requied for sampling larger data sets e.g. when nsam = 1000
  time<-system.time(
    {
      if(!require(PopGenome)){
        install.packages("PopGenome")
      }
      pathtodata<-paste("../../data/",repo,"/",sep="")
      thetas <- seq(theta_lower, theta_upper, increment)
      df <- data.frame()
      initial <- TRUE
      for(t in thetas){
        print(paste("iteration: theta = ",t,sep = ""))
        # reading in data
        data <- NULL
        # print('data is ')
        # set.seed(101)
        # print(data)
        data <- readMS(paste(pathtodata,nsam,"_",nrep,"_",t,".txt",sep=""))
        
        # print(data@n.sites)
        # compute statistics 
        data<-neutrality.stats(data)
        data <- F_ST.stats(data)
        data <- detail.stats(data)
        
        # query statistics
        theta <- rep(t,nrep)
        S_n <- data@n.segregating.sites[,1]
        H_n <- c()
        for(i in 1:nrep){
          H_n <- c(H_n, length(data@region.stats@haplotype.counts[[i]]))
        }
        Pi_n <- data@nuc.diversity.within[,1]
        H_diversity <- c()
        for(i in 1:nrep){
          H_diversity <- c(H_diversity,data@region.stats@haplotype.diversity[[i]][,1][[1]])
        }
        # this is filled in O(mn)
        Eta <- matrix(0,nrep,ceiling((nsam/2))) # a matrix
        for(i in 1:nrep){
          minor <- round(data@region.stats@minor.allele.freqs[[i]]*nsam)
          for(j in 1:length(minor)){
            Eta[i,minor[j]] <- Eta[i,minor[j]] + 1 
          }
        }
        
        # neutrality test statistics can be recorded but wth, they are gonna mean at 0 anyway 
        
        
        # storing the data 
        new_data <- data.frame(theta,S_n,H_n,Pi_n,H_diversity,Eta)
        # write.csv(new_data,paste('temp',t,'.csv',sep=''),row.names=FALSE)
        # print(new_data)
        if(initial){
          df <- new_data
          initial <- FALSE
        } else {
          df <- rbind(df,new_data)
        }
        
      }
    }
  )
  
  write.csv(df,paste(pathtodata,nsam,"_",nrep,".csv",sep=""),row.names=FALSE)
  
  return(time)
}

# let get this on
# df <- compute_ss(100,100,10,150,"repo2")

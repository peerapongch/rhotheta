# data simulation 
# from none
# to summary_stats.R

simulate_ms <- function(nsam, nrep, theta_lower, theta_upper, increment, repo){
  # setwd("C:/cs350/msdir")
  thetas <- seq(theta_lower, theta_upper, increment)
  data_path <- paste("../../data/",repo,"/",sep="")
  for(t in thetas){
    filename <- paste(data_path,nsam,"_",nrep,"_",t,".txt",sep="")
    cmd <- paste("bash -c \" ../../msdir/ms ", nsam, " ", nrep, ' -seeds ', t, ' ', t, " ", t, " -t ", t, " > ", filename, "\"", sep = "")
    print(cmd)
    system(cmd)
  }
}

simulate_ms_rho<-function(nsam, nrep, theta_lower, theta_upper, increment, repo, rho, nsites){
	# setwd("C:/cs350/msdir")
	thetas <- seq(theta_lower, theta_upper, increment)
  data_path <- paste("../../data/",repo,"/",sep="")
  for(t in thetas){
    filename <- paste(data_path,nsam,"_",nrep,"_",t,".txt",sep="")
    cmd <- paste("bash -c \" ../../msdir/ms ", nsam, " ", nrep, ' -seeds ', t, ' ', t, " ", t, " -t ", t, 
    	' -r ', rho, ' ',nsites,
    	" > ", filename, "\"", sep = "")
    print(cmd)
    system(cmd)
  }
}
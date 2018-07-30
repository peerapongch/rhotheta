# simulate under neutral model with mutation and recombination 
nsam<-197
nrep<-100
theta_lower<-10
theta_upper<-200
increment<-10
rho<-387
suffix<-paste('rho',rho,sep='')
repo<-paste('training_',suffix,sep='')
nsites<-5000 #equivalent to windows size

source('data_simulate_ms.R')
simulate_ms_rho(nsam, nrep, theta_lower, theta_upper, increment, repo, rho, nsites)

system.time(
  for (i in c(1)){
    source('data_summary_stats.R')
    compute_ss(nsam, nrep, theta_lower, theta_upper, increment, repo)
  }
)

source('data_train_test_split.R')
data_train_test_split(nsam,nrep,repo,0,suffix)

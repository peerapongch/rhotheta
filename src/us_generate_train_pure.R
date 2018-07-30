# simulate under nuetral model with theta alone 
nsam<-1000
nrep<-100
theta_lower<-10
theta_upper<-150
increment<-10
repo<-'training_pure'
suffix<-'pure'

source('data_simulate_ms.R')
simulate_ms(nsam, nrep, theta_lower, theta_upper, increment, repo)

system.time(
  for (i in c(1)){
    source('data_summary_stats.R')
    compute_ss(nsam, nrep, theta_lower, theta_upper, increment, repo)
  }
)

source('data_train_test_split.R')
data_train_test_split(nsam,nrep,repo,0,suffix)

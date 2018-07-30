# note that the script is written for mutation rate per nwindow base pairs
library(PopGenome)
# for loop over all chromosomes, remember to change all the frompos and topos
pop<-'ZI'
tid<-'2L'
frompos<-'7500000'
topos<-'7600000'
foldername<-paste('../../data/genomes/',pop,'_Chr',tid,'_',frompos,'_',topos,'_fasta',sep='')
# file<-paste(foldername,'/',pop,'_Chr',tid,'_',frompos,'_',topos,'.fasta',sep='')
width<-'5000'
jump<-'5000'
nwindows<-ceiling((as.numeric(topos)-as.numeric(frompos))/as.numeric(width))

system.time(
  for(i in c(1)){
    data<-readData(foldername,format='fasta',include.unknown = FALSE, # this is a very misleading option
                   big.data=TRUE,FAST=TRUE) 
    data.windows<-sliding.window.transform(data,width=as.numeric(width),jump=as.numeric(jump),type=2)
    # data<-read.big.fasta(file,window=as.numeric(width),big.data=TRUE,FAST=TRUE)
  }
)

nsam<-length(get.individuals(data)[[1]])
# non-overlapping sliding windows analysis
## compute
data.windows<-count.unknowns(data.windows)
dfw<-neutrality.stats(data.windows)
dfw<-diversity.stats(dfw,pi=TRUE)
dfw<-F_ST.stats(dfw)
dfw<-detail.stats(dfw)

## query stats for boosting models
S_n<-as.data.frame(dfw@n.segregating.sites)[[1]]
H_n <- c()
for(i in 1:nwindows){
	# reference sequence excluded
	H_n <- c(H_n, length(dfw@region.stats@haplotype.counts[[i]]))
}
Pi_n <- as.data.frame(dfw@nuc.diversity.within)[[1]]
H_diversity <- c()
for(i in 1:nwindows){
  value<-dfw@region.stats@haplotype.diversity[[i]][,1][[1]] # what a bug
  if(is.null(value)){
    value<-0
  }
  H_diversity <- c(H_diversity, value)
}
Eta <- matrix(0,nwindows,ceiling((nsam)/2)) # a matrix
for(i in 1:nwindows){
  minor <- round(dfw@region.stats@minor.allele.freqs[[i]]*nsam) # fix rounding error
  for(j in 1:length(minor)){
    Eta[i,minor[j]] <- Eta[i,minor[j]] + 1 
  }
}

# query analytic estimation statistics
theta_watterson<-as.data.frame(dfw@theta_Watterson)[[1]]
theta_tajima<-as.data.frame(dfw@theta_Tajima)[[1]]
theta_fu<-as.data.frame(dfw@theta_Fu.Li)[[1]] # need to confirm what this value corresponds to

# compute theta_eta1
mult<-(nsam-1)/nsam
theta_eta1<-Eta[,1]*mult

# dont forget the number of sites 
bp<-round((1-dfw@missing.freqs)*as.numeric(width))[,1]
# region<-row.names(dfw@theta_Watterson)

# adding some test statistics 
test_tajima_d<-dfw@Tajima.D[,1]
test_fu_li_d<-dfw@Fu.Li.D[,1]
test_fu_li_f<-dfw@Fu.Li.F[,1]

# wrap up and save
# omit region
new_data <- data.frame(S_n,H_n,Pi_n,H_diversity,Eta,
                       theta_watterson,theta_tajima,theta_eta1,theta_fu,bp,
                       test_tajima_d,test_fu_li_d,test_fu_li_f)
write.csv(new_data,paste('../../data/genomes/',pop,'_Chr',tid,'_',frompos,'_',topos,'_window',width,'.csv',sep=''),row.names=FALSE)

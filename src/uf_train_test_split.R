# this script standardises the input data then split them into train and test X and y sets before saving variables to rda 
# from summary_stats.R
# to various


data_train_test_split<-function(nsam,nrep,repo,nattributes,suffix){
	require(caret)
	data_path <- paste("../../data/",repo,"/",sep="")
	data <- read.csv(paste(data_path,nsam,"_",nrep,".csv",sep=""))
	# if (nattributes!=0){
	# 	data <- data[,1:(nattributes+1)] # customisation valid after boosting found significant variables
	# }
	data<-data[,colSums(data)!=0]
	# partitioning data
	inTraining <- createDataPartition(data$theta, p=0.70, list=FALSE)
	train <- data[inTraining,]
	test <- data[-inTraining,]
	rownames(test) <- seq(length=nrow(test))
	rownames(train) <- seq(length=nrow(train))

	# keep mean and sd of each column
	mus_train <- apply(X=train[,2:dim(train)[2]] ,FUN=mean,MARGIN=2)
	sigmas_train <- apply(X=train[,2:dim(train)[2]] ,FUN=sd,MARGIN=2)
	# mus_test <- apply(X=test[,2:dim(test)[2]],FUN=mean,MARGIN=2)
	# sigmas_test <- apply(X=test[,2:dim(test)[2]],FUN=sd,MARGIN=2)

	# centre and scale data 
	train[,2:dim(train)[2]] <- scale(train[,2:dim(train)[2]], center=TRUE, scale=TRUE)
	test[,2:dim(test)[2]] <- scale(test[,2:dim(test)[2]], center=mus_train, scale=sigmas_train)
	# # 'fixing NAs' by essentially drop them ==> this is not good
	# train <- train[,colSums(is.na(train))<nrow(train)]
	# test <- test[,colSums(is.na(test))<nrow(test)]
	## alternatively
	train[,colSums(is.na(train))>=nrow(train)]<-0
	test[,colSums(is.na(test))>=nrow(test)]<-0
	# separate covariates and response 
	train.X <- train[,2:dim(train)[2]]
	train.Y <- train[,1]
	test.X <- test[,2:dim(test)[2]]
	test.Y <- test[,1]

	# additional for training boosted gam and tree 
	columns<-colnames(train[,colSums(train)!=0])
	columns.X<-colnames(train.X[,colSums(train.X)!=0])

	save(test.Y,test.X,train.X,train.Y,train,test,mus_train,sigmas_train,columns,columns.X,
		file=paste('../../data/rdas/traindata_',nsam,'_',nrep,'_',suffix,'.rda',sep=''))
}


testdata_scale<-function(nsam,nrep,repo,train_data){
	# print(getwd())
	load(train_data) # this better points to rdas folder
	data_path <- paste("../../data/",repo,"/",sep="")
	test <- read.csv(paste(data_path,nsam,"_",nrep,".csv",sep=""))
	test[,2:dim(test)[2]]<-t(apply(X=test[,2:dim(test)[2]],FUN=function(x){(x-mus_train)/sigmas_train},MARGIN=1))
	# test[,2:dim(test)[2]]<-t(apply(X=test[,2:dim(test)[2]],FUN=function(x){x/sigmas_train},MARGIN=1))
	test.X <- test[,2:dim(test)[2]]
	test.Y <- test[,1]
	save(test,test.X,test.Y,file=paste("../../data/rdas/testdata_",nsam,'_',nrep,'.rda',sep=''))
}

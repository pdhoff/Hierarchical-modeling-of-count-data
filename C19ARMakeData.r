## -- weekly county-level data on numbers of reported cases
#Ccases<-readRDS(url("https://github.com/pdhoff/US-counties-C19-data/blob/master/UScountiesC19Cases.rds?raw=true"))
Ccases<-readRDS("~/Dropbox/GHDBlog/COVID19-Data/UScountiesC19Cases.rds")


## -- weekly county-level data on number of reported deaths
#Cdeaths<-readRDS(url("https://github.com/pdhoff/US-counties-C19-data/blob/master/UScountiesC19Deaths.rds?raw=true")) 
Cdeaths<-readRDS("~/Dropbox/GHDBlog/COVID19-Data/UScountiesC19Deaths.rds")


## -- convert to state level data
stateFIPS<-substring(rownames(Ccases),1,2)

Scases<-Sdeaths<-NULL
for(st in sort(unique(stateFIPS))){

  Scases<-rbind(Scases, apply(Ccases[ stateFIPS==st,,drop=FALSE],2,sum ) )
  Sdeaths<-rbind(Sdeaths, apply(Cdeaths[ stateFIPS==st,,drop=FALSE],2,sum ) )  
}
rownames(Scases)<-rownames(Sdeaths)<-sort(unique(stateFIPS)) 
colnames(Scases)<-colnames(Sdeaths)<-colnames(Ccases)


## -- state populations 
USCdata<-readRDS(url("https://github.com/pdhoff/US-counties-data/blob/master/UScounties.rds?raw=true"))
Spops<-tapply(USCdata$population, USCdata$stateFips,sum)
all(names(Spops) == rownames(Scases))


## --  design matrix for an AR fit: 
yLaggedyx<-function(y,x,lag=1,loglag=TRUE){

  # y is the outcome time series that will be lagged 
  # x is a predictor time series that will be lagged

  LX<-LY<-NULL
  for(l in 1:lag){
    LX<-cbind(LX, c(rep(NA,l),x[1:(length(x)-l)] )  )
    LY<-cbind(LY, c(rep(NA,l),y[1:(length(y)-l)] )  )
  }
  if(loglag){ LX<-log(LX+1) ; LY<-log(LY+1) }  
  yX<-cbind(y,LY,LX) 
  yX<-yX[!apply(is.na(yX),1,any ) ,]
  list(y=yX[,1],X=yX[,-1] )
}

## -- choose lags for modeling
lag<-4


## -- construct design matrices   
yd<-yc<-X<-W<-list()
for(i in 1:nrow(Sdeaths)){ 
  yX<-yLaggedyx(Sdeaths[i,],Scases[i,],lag=lag)   
  yd[[i]]<-yX$y ; X[[i]]<-yX$X 
  colnames(X[[i]])<-c(paste0( "dlag.m",1:lag),paste0( "clag.m",1:lag))
  yc[[i]]<-tail(Scases[i,],length(yd[[1]])) 
   
  ## predictors for coefficients
  Wi<-cbind(0,diag(nrow=1+ncol(X[[1]])))  
  Wi[1,1:2]<-c(1,log(Spops[i]))  
  W[[i]]<-Wi
}
names(yd)<-names(yc)<-names(X)<-names(W)<-unique(USCdata$state)

save(yd,yc,X,W,file=paste0("C19AR",lag,".rdata") )



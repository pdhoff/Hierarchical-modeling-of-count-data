## -- functions for downloading and wrangling data 
source("https://raw.githubusercontent.com/pdhoff/US-counties-C19-data/master/USC19data.r") 


## -- download and convert data to weekly state-level counts 
X<-stateify( weekify( pullC19data() ) )
Sdeaths<-X[,,1] 
Scases<-X[,,2]


## -- state populations 
USCdata<-readRDS(url("https://github.com/pdhoff/US-counties-data/blob/master/UScounties.rds?raw=true"))
Spops<-tapply(USCdata$population, USCdata$stateFips,sum)

## -- check 
all( fips2state(names(Spops))  == rownames(Scases))


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



---
title: "State-level C19 analyses"
author: "Peter Hoff"
date: "June 17, 2020"
output: github_document
---

### Summary 
This document illustrates the use of a hierarchical negative binomial model 
to analyze time-series of US state-level C19 counts of deaths and reported cases. 
This is done using essentially a bivariate fourth-order autoregressive generalized linear model. For example, the death count $y^d_{j,t}$ in state $t$ and week $j$ is modeled as a function of the (log) death counts and case counts in the previous four weeks: 
\[
 y^d_{j,t}  \sim   \beta_{j,0} + \beta_{j,1} \log y^d_{j,t-1} + \cdots +
       \beta_{j,4} \log y^d_{j,t-4}  + 
         \beta_{j,5} \log y^c_{j,t-1} + \cdots +
       \beta_{j,8} \log y^c_{j,t-4} .
\]
A similar model is fit for the case data. 



```{r}
## -- load data that has been set up for AR model fitting
lag<-4
load(paste0("C19AR",lag,".rdata")) 
```

```{r}
## -- load functions for hierarchical negative binomial model fitting
source("negbinHGLM.r") 
```

```{r}
## -- fit AR model to deaths and to cases
fitDeaths<-negbinHGLM(yd,X,W)
fitCases<-negbinHGLM(yc,X,W)
```

```{r}
## -- get fitted values for deaths
Dobs<-Dfit<-NULL
for(i in 1:length(yd)){

   Bi<-fitDeaths$BPS[i,,]
   Li<-fitDeaths$LPS[,i]
   Xi<-cbind(1,X[[i]])
   Dfit<-rbind(Dfit,apply(sweep(exp(Xi%*%Bi),2,Li,"*"),1,mean))
   Dobs<-rbind(Dobs,yd[[i]] )
}
```

## -- get fitted values for cases
Cobs<-Cfit<-NULL
for(i in 1:length(yc)){

   Bi<-fitCases$BPS[i,,]
   Li<-fitCases$LPS[,i]
   Xi<-cbind(1,X[[i]])
   Cfit<-rbind(Cfit,apply(sweep(exp(Xi%*%Bi),2,Li,"*"),1,mean))
   Cobs<-rbind(Cobs,yc[[i]] )
}
```

```{r}
## -- plot fitted and observed
par(mfrow=c(5,10),mar=c(0,0,0,0),mgp=c(0,0,0))
for(i in (1:51)[-9]){
 plot(c(1,ncol(Dobs)),c(0,1.1),type="n",xlab="",ylab="",xaxt="n",yaxt="n")

 yo<-Cobs[i,] ; yf<-Cfit[i,] ; mx<-max(c(yo,yf)) 
 yo<-yo/mx ; yf<-yf/mx  
 lines(yf,col="lightblue",lwd=3)
 lines(yo,col="blue") 

 yo<-Dobs[i,] ; yf<-Dfit[i,] ; mx<-max(c(yo,yf)) 
 yo<-yo/mx ; yf<-yf/mx 
 lines(yf,col="pink",lwd=3)
 lines(yo,col="red")  
 text(1,max(c(yo,yf))*1.05,names(X)[i],pos=4,cex=.8 ) 
}
```

## -- forecasting 

# - function to generate forcasts from the two time series models
nbforecast<-function(deaths,cases,dbeta,cbeta,dlambda,clambda,nf=4){

  lag<-(length(dbeta)-1)/2

  for(i in 1:nf)
  {
    xf<-c(1,log( c( tail( deaths,lag), tail(cases,lag) ) +1 ))
    lpd<-(xf%*%dbeta)
    lpc<-(xf%*%cbeta)
    
    d<-rnbinom(1,mu=dlambda*exp(lpd),size=dlambda)
    c<-rnbinom(1,mu=clambda*exp(lpc),size=clambda)
    
    deaths<-c(deaths,d)
    cases<-c(cases,c)
  }
  cbind( tail(deaths,lag), tail(cases,lag)  )
}

 
## -- Now simulate from the forecast distribution 5 times for 
## -- each set of parameter values simulated from the posterior 
## -- distribution. This accounts for uncertainty in the model 
## -- parametre estimates. 

nf<-4 
nsim<-nrow(fitDeaths$LPS)*5
sim<-rep(1:nrow(fitDeaths$LPS),5)  

FCAST<-array(dim=c(length(yd),nf,nsim,2) ) 
for(i in 1:length(yd))
{ 
  dforecast<-NULL 
  for(s in 1:nsim)
  {  
   ss<-sim[s] 
   dbeta<-fitDeaths$BPS[i,,ss] 
   cbeta<-fitCases$BPS[i,,ss] 
   dlambda<-fitDeaths$LPS[ss,i] 
   clambda<-fitCases$LPS[ss,i] 
   fi<-nbforecast(yd[[i]],yc[[i]],dbeta,cbeta,dlambda,clambda,nf)
   FCAST[i,,s,]<-fi  
  }
}



## -- Compute point (median) forecast and glom onto fitted values
MFCAST<-apply(FCAST,c(1,2,4),median)  
DFF<-cbind(Dfit,MFCAST[,,1] ) 
CFF<-cbind(Cfit,MFCAST[,,2] )

## -- Plot fitted, observed and forecasted values all together
par(mfrow=c(5,10),mar=c(0,0,0,0),mgp=c(0,0,0))
for(i in (1:51)[-9]){
 plot(c(1,ncol(DFF)),c(0,1.1),type="n",xlab="",ylab="",xaxt="n",yaxt="n")

 yo<-Cobs[i,] ; yf<-CFF[i,] ; mx<-max(c(yo,yf))
 yo<-yo/mx ; yf<-yf/mx
 lines(yf,col="lightblue",lwd=3)
 lines(yo,col="blue")

 yo<-Dobs[i,] ; yf<-DFF[i,] ; mx<-max(c(yo,yf))
 yo<-yo/mx ; yf<-yf/mx
 lines(yf,col="pink",lwd=3)
 lines(yo,col="red")
 text(1,max(c(yo,yf))*1.05,names(X)[i],pos=4,cex=.8 )
}



rbetaGibbs<-function(y,X,beta,lambda,
            mu=rep(0,length(beta)),iSigma=diag(nrow=length(beta)))  
{ 
  # Simulation of regression coefficients of a neg binom regression model

  # y is the vector of negative binomial outcomes 
  # X is the design matrix
  # beta is the current value of beta
  # lambda is the size parameter 
  # mu is the expectation of beta  
  # iSigma is the inverse-variance of beta
  
  XW<-sweep(X,1,BayesLogit::rpg(length(y),lambda+y, X%*%beta),"*")  
  Q<-crossprod(X,XW)
  l<-crossprod(X,(y-lambda)/2 ) 
  V<-solve( Q + iSigma )
  m<-V%*%( iSigma%*%mu + l )
  m + t(chol(V))%*%rnorm(length(l)) 
}


rlambdaMH<-function(y,X,beta,lambda,a,b)
{
  # MH-update of size parameter in negative binomial regression

  # y is the vector of negative binomial outcomes 
  # X is the design matrix
  # beta is the current value of beta
  # lambda is the size parameter
  # (a,b) are parameters of the inverse-gamma prior for lambda

  eXb<-exp(X%*%beta)
  lambda0<-lambda
  mu0<-lambda0*eXb 

  lambda1<-1/exp( log(1/lambda) + rnorm(1,0,.25) ) 
  mu1<-lambda1*eXb

  lhr<-sum( dnbinom(y,size=lambda1,mu=mu1,log=TRUE)   -  
            dnbinom(y,size=lambda0,mu=mu0,log=TRUE) ) + 
            dgamma(1/lambda1,a,b,log=TRUE) - dgamma(1/lambda,a,b,log=TRUE) + 
            log(1/lambda1) - log(1/lambda)

  if(log(runif(1))<lhr){ lambda<-lambda1 }
  lambda
}


rlambdaGA<-function(y,X,beta,lambda,a,b)
{ 
  # Grid update of size parameter in negative binomial regression

  # y is the vector of negative binomial outcomes 
  # X is the design matrix
  # beta is the current value of beta
  # lambda is the size parameter
  # (a,b) are parameters of the inverse-gamma prior for lambda

  lambdas<-seq(.5*lambda,1.5*lambda,length=100) 
  lgpy<-apply(lgamma( outer(lambdas,y,"+") ),1,sum) 
  lpl<- lgpy - length(y)*lgamma(lambdas) - lambdas*sum(log(1+exp(X%*%beta))) +
        -(a+1)*log(lambdas) - b/lambdas  
  sample(lambdas,1,prob=exp(lpl-max(lpl))) 
}



negbinHGLM<-function(y,X,W=NULL,nmc=5000,nburn=1000,thin=10,plot=TRUE,seed=1)
{

  # MCMC for a hierarchical negative binomial model of the form
  #
  # y[[i]] ~ negative binomial( size=Lambda[i],mu=Lambda[i]*exp(X[[i]]%*%B[i,])
  # B[i,] ~ N( W[[i]]%*%alpha , Sigma ) 
  # 1/Lambda[i] ~ gamma(a,b) 
 
  # y, X and W are lists of group-specific outcomes, design matrices, 
  # and design matrices for the regression parameters.
  # So the number of rows of X should match the length 
  # of y, and the number of rows of W should be the number of 
  # columns of X. An intercept will be added to X if it is not 
  # already there. 

  # nmc - the number of MCMC iterations, after burn in. 
  # nburn - the number of burn in iterations.
  # thin - the frequency with which simulated parameters are saved.  
  # plot - indicator if results should be plotted during the MCMC iterations 
  # seed - a random seed


 

  ## random seed
  set.seed(seed)


  ## add intercept to X if needed
  if( ! all( sapply(X,function(x){ all(x[,1]==1) }) )    )
  {
    for(i in 1:length(X)){ X[[i]]<-cbind(1,X[[i]]) }  
  } 


  ## construct W matrices if needed 
  if(is.null(W))
  { 
    W<-list() 
    for(i in 1:length(y)){ W[[i]]<-diag(nrow=ncol(X[[1]])) }
  }


  ## starting values 
  Beta<-NULL ; Lambda<-NULL ; Fail<-NULL
  for(i in 1:length(y))
  {
    fit<-suppressWarnings(MASS::glm.nb(y[[i]]~X[[i]][,-1]))
    Fail<-c(Fail, any(zapsmall(fit$fitted.values)==0)) 
    beta<-fit$coef ; lambda<-fit$theta 
    beta[1]<-beta[1]-log(fit$theta)
    Beta<-rbind(Beta,beta)
    Lambda<-c(Lambda,fit$theta)  
    colnames(Beta)<-colnames(X[[1]]) 
  }
  if(sum(Fail)>0)
  { 
    Lambda[Fail]<-median(Lambda[!Fail]) 
    Beta[Fail,]<-0
    Beta[Fail,1]<- log( sapply(y[Fail],mean)/Lambda[Fail] )
  }

  ## starting values for hierarchical params
  alpha<-rep(0,ncol(W[[1]])) ; names(alpha)<-colnames(W[[1]]) 
  iSigma<-diag(nrow=ncol(Beta))  

  ## empirical Bayes estimates of inverse gamma distribution for size params
  mllab<-function(ab){ -sum(dgamma(1/Lambda,exp(ab[1]),exp(ab[2]),log=TRUE)) }
  ab<-exp(optim(c(1,1),mllab)$par) 
  a<-ab[1] ; b<-ab[2]


  ## -- MCMC  

  # - objects where parameter values will be saved
  BPS<-array(dim=c(length(y),ncol(Beta),nmc/thin)) ; LPS<-NULL 
  APS<-NULL ; iSPS<-matrix(0,ncol(Beta),ncol(Beta)) ; ABPS<-NULL 

  # - progress bar 
  progressBar<-txtProgressBar(0,1)

 
  # - MCMC loop 
  for(s in -nburn:nmc)
  {  

    ## -- update group level parameters  
    for(i in 1:length(y))
    {
      Beta[i,]<-rbetaGibbs(y[[i]],X[[i]],Beta[i,],Lambda[i],
                           mu=W[[i]]%*%alpha,iSigma=iSigma)
      Lambda[i]<-rlambdaGA(y[[i]],X[[i]],Beta[i,],Lambda[i],a,b)
    }

    ## -- update hierarchical model for regression coef  

    # - update alpha
    Q<-diag(nrow=length(alpha))/1000 
    l<-rep(0,length(alpha)) 
    for(i in 1:length(y))
    {  
      Q<-Q+t(W[[i]])%*%iSigma%*%W[[i]]
      l<-l+t(W[[i]])%*%iSigma%*%Beta[i,] 
    }
    Valpha<-solve(Q) 
    Ealpha<-Valpha%*%l 
    alpha<- Ealpha+t(chol(Valpha))%*%rnorm(length(alpha))  

    # - update Sigma
    SS<-diag(ncol(Beta))
    for(i in 1:length(y)){ SS<-SS+tcrossprod( Beta[i,]- W[[i]]%*%alpha ) }
    iSigma<-rWishart(1,ncol(Beta)+2+nrow(Beta), solve(SS) )[,,1]


    ## -- output
    if(s%%thin==0)
    {  
      setTxtProgressBar(progressBar,(s+nburn)/(nburn+nmc+1) ) 

      if(s>0)
      {
        BPS[,,s/thin]<-Beta 
        LPS<-rbind(LPS,Lambda) 
        APS<-rbind(APS,c(alpha))
        iSPS<-iSPS+iSigma 

        if(plot)
        { 
          par(mfrow=c(1,3),mar=c(3,3,1,1),mgp=c(1.75,.75,0)) 
          boxplot(Beta,ylab=expression(beta)) ;abline(h=0,col="gray",lty=2)  
          matplot(APS,ylab=expression(alpha),type="l",lty=1)    
          hist(1/Lambda,prob=TRUE,main="",ylab="",xlab=expression(1/lambda)) 
          x<-seq(min(1/Lambda),max(1/Lambda),length=100)  
          lines(x,dgamma(x,a,b)) 
        }
      }
      close(progressBar) 
    }
  }
  list(BPS=BPS,LPS=LPS,APS=APS,SigmaHat=solve(iSPS/nrow(LPS)),ab=ab)
}



library("data.table")
## library("Rcmdr")
library(ggplot2)
library(reshape2)
library("DEoptim")

#######
## Wrapper functions
#######

## Function to optimize model
optim.mod <- function(response,design.mat,ll.fun,parms,lower=NULL,upper=NULL,restarts=200,
                      parm.min=-10,parm.max=10,method=NULL,error.mod,maxit=NULL,refine=FALSE){
  ## make sure inputs are in matrix form
  response=as.matrix(response)
  design.mat=as.matrix(design.mat)
  ## Initialize parameters
  nparms=length(parms)
  optim.fun="optim"
  opt=list(par=parms,fn = ll.fun,response=response,design.mat=design.mat,
           error.mod=error.mod)
  if (!is.null(method)){
    opt[["method"]]=method
  }  else if(!is.null(lower) || !is.null(upper)){
    opt[["method"]]="L-BFGS-B"
    if(!is.null(lower)){
      opt[["lower"]]=lower  
    } 
    if(!is.null(upper)){
      opt[["upper"]]=upper     
    }
  } else {
    opt[["method"]]="Nelder-Mead"
  }
  
  if(!is.null(maxit)){
    opt[["control"]]=list(maxit=maxit)
  }
  
  
  ## Run the optimizer of choice
  best=do.call(optim.fun,args = opt)
  if(restarts>1){
    for(j in 2:restarts){
      opt[["par"]]=runif(n = nparms,parm.min,parm.max)    
      mod=do.call(optim.fun,args = opt)
      if(mod$value<best$value){
        best=mod
      }  
    }
  }
  ## If refine, run one additional gradient descent
  if(refine){
    opt.ref=list(par=best$par,fn = ll.fun,response=response,design.mat=design.mat,
                 error.mod=error.mod,control=list(maxit=3000))
    best=do.call(optim,args = opt.ref)
  }
  return(best) 
}

## Optimization using evolutionary algorithm
deOptim.mod<-function(response,design.mat,ll.fun,parms,lower,upper,restarts=10,
              error.mod,refine=FALSE){
  response=as.matrix(response)
  design.mat=as.matrix(design.mat)
  ## Initialize parameters
  nparms=length(parms)
  if(length(lower)==1){
    lower=rep(lower,nparms)
  } 
  if(length(upper)==1){
    upper=rep(upper,nparms)
  }
  opt=list(fn = ll.fun,lower=lower,upper=upper,response=response,design.mat=design.mat,
           error.mod=error.mod,control=list(itermax=10000,NP=20*length(lower),trace=100))
  best=do.call(DEoptim,args = opt)
  if(refine){
      best=optim.mod(response = response,design.mat = design.mat, ll.fun=ll.fun, 
                parms=as.numeric(best$optim$bestmem),error.mod=error.mod,
                restarts=1,refine=FALSE)
    }
  return(best)
}
 
## Take model and get predicted values for specific inputs
model.pred <- function(design.mat, params, val.fun){
  design.mat=as.matrix(design.mat)
  do.call(val.fun,args = list(design.mat=design.mat,x=params))  
}

####
## Specific models
####

##
## Linear models
##

## Linear act
linear.act <- function(x,design.mat){
  act=x[1] + as.matrix(design.mat) %*% as.matrix(x[2:(1+ncol(design.mat))])
  return(act)
}

## Linear model
linear.mod <- function(x,response,design.mat,
                       error.mod=c("gaussian","log-normal")){
  if(length(error.mod)>1){
    stop("Must choose an error model.")
  }
  act=linear.act(x,design.mat)
  expression=linear.val(act)
  if(error.mod=="gaussian"){
    ll=dnorm(response,act,x[2+ncol(design.mat)],log=TRUE)    
  } else if(error.mod=="log-normal"){
    ll=dlnorm(response,my.log(expression),x[2+ncol(design.mat)],log=TRUE)
  } else {
    stop("Not a valid error model.")
  }
  sum.ll=sum(ll)
  return(-sum.ll)
}

## Linear value 
linear.val <- function(act,x=NULL){
  return(act)
}

## Returns values for a set
linear.predict <- function(design.mat,x){
  act=linear.act(x,design.mat)
  express=linear.val(act)
  return(data.table(design.mat,activity=as.numeric(act),expression=as.numeric(express)))  
}

##
## Multiplicative models
##


multiplicative.model <- function(x,response,design.mat,
                                 error.mod=c("gaussian","log-normal","gamma")){
  act=linear.act(x,design.mat)
  expression=multiplicative.val(act)
  if(error.mod=="gaussian"){
    ll=dnorm(response,expression,x[2+ncol(design.mat)],log=TRUE)
  } else if(error.mod=="log-normal"){
    ll=dlnorm(response,my.log(expression),x[2+ncol(design.mat)],log = TRUE)
  } else if(error.mod=="gamma"){
    ll=dgamma(response,expression/x[2+ncol(design.mat)],scale = x[2+ncol(design.mat)],
              log=TRUE)
  } else {
    stop("Not a valid error model.")
  }
  sum.ll=sum(ll)
  return(-sum.ll)
}
## Set up in independant (multiplicative) model
multiplicative.val <- function(act,x=NULL){
  return(exp(act))
}

## Get predicted model value given parameters and design matrix
multiplicative.predict <- function(design.mat,x){
  act=linear.act(x,design.mat)
  expression=multiplicative.val(act)
  return(data.table(design.mat,activity=as.numeric(act),
                    expression=as.numeric(expression)))  
}

##
## Logistic models
##

logistic.model <- function(x,response,design.mat,
                           error.mod=c("gaussian","scaled-gaussian",
                                       "log-normal","gamma")){
  act=linear.act(x,design.mat)
  expression=logistic.val(act,x[ncol(design.mat)+2])
  if(error.mod=="gaussian"){
    ll=dnorm(response,expression,x[3+ncol(design.mat)],log=TRUE)
  } else if(error.mod=="log-normal"){
    ll=dlnorm(response,my.log(expression),x[3+ncol(design.mat)],log = TRUE)
  } else if(error.mod=="gamma"){
    ll=dgamma(response,expression/x[3+ncol(design.mat)],scale = x[3+ncol(design.mat)],
              log=TRUE)
  } else if(error.mod=="scaled-gaussian"){
    ll=unlist(apply(cbind(response,expression),1,function(y) 
      dnorm(y[1],y[2],y[2]*x[3+ncol(design.mat)],log=TRUE)))
  }else {
    stop("Not a valid error model.")
  }
  sum.ll=sum(ll)
  return(-sum.ll)
}

logistic.val <- function(act,scale){
  return(scale/(1+exp(-act)))
}

## Get predicted model value given parameters and design matrix
logistic.predict <- function(design.mat,x){
  act=linear.act(x,design.mat)
  expression=logistic.val(act,x[ncol(design.mat)+2])
  return(data.table(design.mat,activity=as.numeric(act),
                    expression=as.numeric(expression)))  
}

#####
## Plotting functions
#####
## The dots are any aditional parameters that must be passed to val.fun
plot.model <- function(model,observed,design.mat,val.fun,error.type,...){
  ## Get activity values for actual oberservations
  act=linear.act(model$par,as.matrix(design.mat))
  real=data.table(activity=as.numeric(act),observed=observed)
  
  sim.act=seq(min(real$activity)-abs(min(real$activity))*0.1,max(real$activity)+abs(max(real$activity))*0.1,by=0.01)
  dots=list(...)
  dots$act=sim.act
  sim.exp=do.call(val.fun,args = dots)
  out=data.frame(activity=sim.act,expression=sim.exp)
  
  if(error.type=="log-normal"){
    err=lognormal.error.int(out$activity,out$expression,sdlog = model$par[length(model$par)],quart = c(0.05))
  } else if (error.type=="gaussian"){
    err=gaussian.error.int(out$activity,out$expression,sd = model$par[length(model$par)],quart = c(0.05))    
  } else{
    stop("Error type not implemented")
  }
  cc <- scales::seq_gradient_pal("light blue", "blue", "Lab")(seq(0,0.5,length.out=length(unique(err$Quantile))))
  
  g=ggplot()+
    geom_polygon(data=err,aes(x=x,y=y,fill=Quantile),alpha=1)+
    scale_fill_manual(values=cc)+
    geom_line(data=out,aes(activity,expression),color="black",size=2)+
    geom_point(data=real,aes(activity,observed),color="red")+
    theme_bw(base_size = 28)+
    xlab("Activity/-Energy")+
    ylab("Expression")+
    theme(legend.position = c(0, 1),legend.justification = c(0, 1))
  
  return(g)
}

lognormal.error.int<-function(activity,expression,sdlog,quart){
  if(sum(quart>0.5)>0){
    stop("Quantiles 0<quart<0.5")
  }
  err=list()
  for(q in quart){
    err[[as.character(q)]]=data.table(x=c(activity,rev(activity)),
                                      y=c(qlnorm(q, log(expression), sdlog),
                                          rev(qlnorm(1-q, log(expression), sdlog))),
                                      Quantile=paste0(q,"-",1-q))
  }
  return(rbindlist(err))
}

gaussian.error.int <- function(activity,expression,sd,quart){
  if(sum(quart>0.5)>0){
    stop("Quantiles 0<quart<0.5")
  }
  err=list()
  for(q in quart){
    err[[as.character(q)]]=data.table(x=c(activity,rev(activity)),
                                      y=c(qnorm(q, expression, sd),
                                          rev(qnorm(1-q, expression, sd))),
                                      Quantile=paste0(q,"-",1-q))
  }
  return(rbindlist(err))
}

## Predicted vs. actual plot 
pvr.plot <- function(predicted,observed,condition){
  pvr.dat=data.table(c=as.factor(condition),p=as.numeric(predicted),o=as.numeric(observed))
  means=pvr.dat[,list(p=mean(p),o=mean(o)),by=c]
  g=ggplot()+
    geom_point(data=means,aes(x=p,y=o),color="red",size=5)+
    geom_point(data=pvr.dat,aes(x=p,y=o))+
    geom_abline(slope=1)+
    ylab("Observed")+
    xlab("Predicted")+
    theme_bw(base_size = 20)     
  return(g)
}

plot.residual <- function(observed, predicted, model=NA, logscale=TRUE){
  if(!logscale){
    dat=data.table(p=predicted,o=observed,r=observed-predicted)
    rs="Residual"
  } else{
    dat=data.table(p=log(predicted),o=log(observed),r=log(observed)-log(predicted))
    rs="log(Residual)"    
  }
  if(is.na(model)){
    g=ggplot(dat,aes(x=p,y=r))+
      geom_point()+
      geom_abline(intercept=0,slope=0)+
      theme_bw(base_size = 12)+
      xlab("Model")+
      ylab(rs)
  } else {
    dat[,m:=model]
    g=ggplot(dat,aes(x=m,y=r,group=m,fill=m))+
      geom_boxplot(notch=TRUE)+
      geom_abline(intercept=0,slope=0)+
      theme_bw(base_size = 20)+
      xlab("Model")+
      ylab(rs)
  }
  return(g)
}

## BIC function
bic <- function(fit.model,N){
  return(-2*-fit.model$value+length(fit.model$par)*log(N))
}

my.log<-function(x){
  x[x<0]=-Inf
  x[x==0]=-300
  x[x>0]=log(x[x>0])
  return(x)   
}

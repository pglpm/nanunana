library('magicaxis')
library('mvtnorm')
library('ellipse')
library('LaplacesDemon')
library('RColorBrewer')
library('metRology')
library('magrittr')


densplot <- function (x,adjust=1,...) {
  density(x,adjust) %>%
      plot(.,...)
}

mu <- 3

draws <- rnorm(2,mu,3)
n <- length(draws)

shape <- 2 + n/2
scale <- 1/(1 + sd(draws)*(n-1)/2)

df <- 2*shape
sd <- 1/sqrt(shape*scale)

PGF <- function(data){
    return(rgamma(1,shape=shape,scale=scale))
}


mydata <- list(data=draws, PGF=PGF, mon.names='y', parm.names=c('tau'),N=n, y=mean(draws))

model <- function(parm,number=1){
    return <- rexp(number, rate=parm)
}


hyperprior <- function(parm,data){
    parm <- interval(parm,1e-6,1e6)
    tau.prior <- dgamma(parm,shape=shape,scale=scale, log=T)
    LL <- sum(dnormp(data$data,mean=mu,prec=parm,log=T))
    LP <- LL+ tau.prior
    return <- list(LP=LP, Dev=-2*LL, Monitor=rnormp(1,mean=mu,prec=parm),yhat=rnormp(1,mean=mu,prec=parm),parm=parm)
}

posterior <- function(parm,data){
    LP <- 1
    return <- list(LP=LP,Dev=1,Monitor=1,yhat=rnormp(1,mean=mu,prec=parm),parm=parm)
}


Sample <- LaplacesDemon(hyperprior, mydata, Initial.Values=c(1), Thinning=1,
                     Iterations=10000, Status=200,
                     Algorithm="AFSS",
                     Specs=list(A=500, B=NULL, m=100, n=0, w=1)
                     )

plot(Sample,BurnIn=500,mydata,PDF=TRUE,Parms=NULL)


postpredictive <- predict(Sample,posterior,mydata,CPUs=2)

densplot(postpredictive$yhat,
    adjust=0.1,
    main='test',
    xlab='tau2')

xx <- seq(-8,8,0.1)
yy <- dt.scaled(xx,df,mean=mu,sd=sd)
lines(x=xx,y=yy,col='green')


plot(postpredictive, Style="Density",PDF=TRUE)

stop()
densplot(Sample$Monitor,
    adjust=1,
    main='test',
    xlab='tau')

xx <- seq(-8,8,0.1)
yy <- dt.scaled(xx,df,mean=mu,sd=sd)
lines(x=xx,y=yy,col='green')


quadr <- IterativeQuadrature(hyperprior, c(5), mydata, Covar=NULL, Iterations=500,
                             sir=F, Algorithm="AGHSG",
                             Specs=list(K=5, Kmax=7, Packages=NULL, Dyn.libs=NULL),
                             Stop.Tolerance=c(1e-5,1e-15), CPUs=1)




stop()

## png('testplot2.png')
## magcon(Fit$Posterior2[9000:10000,1],Fit$Posterior2[9000:10000,2],
##        conlevels=c(0.5,0.68,0.95), lty=c(2,1,3), imcol=brewer.pal(n=9,name='Blues'))

## title(xlab='mu',ylab='sigma')

## polygon(ellipse(cov.rob(Fit$Posterior2)$cov,centre=Fit$Summary2[,'Mean'],
##                 level=c(pnorm(1)-pnorm(-1))), col=hsv(v=0,alpha=0.2),border=NA)

## points(mymu,mysigma,col='blue',pch=4)

## points(Fit$Summary2[1,'Median'],Fit$Summary2[2,'Median'],
## col='black',pch=4)
## dev.off()

## png('testplottrajectory2.png')
## magplot(Fit$Posterior1,type='l',col=hsv(alpha=0.3),xlab='mu',ylab='sigma')
## dev.off()


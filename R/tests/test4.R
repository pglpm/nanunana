library('magicaxis')
library('mvtnorm')
library('ellipse')
library('LaplacesDemon')
library('RColorBrewer')

mymu <- 3
mysigma <- 0.1
myrate <- 20

draws <- rnorm(10,-3,mysigma)
n <- length(draws)
print(mean(draws))

print(sd(draws)*(n-1)/n)

mydata <- list(y=c(10), mon.names='y', parm.names=c('theta'),N=1)

model <- function(parm){
    rx <- sample(0:1,1)
    return <- rnorm(1,rx*parm+(1-rx)*(parm+5),mysigma)
}


hyperprior <- function(parm,data){
    parm <- interval(parm,1e-6,1e6)
    theta.prior <- dexp(parm, rate=myrate, log=TRUE)
    rx <- sample(0:1,1)
    LL <- log(dnorm(data$y, parm, mysigma)*0.5 + dnorm(data$y, parm+5, mysigma)*0.5)
    LP <- LL + theta.prior
    return <- list(LP=LP, Dev=-2*LL, Monitor=model(parm),yhat=1,parm=parm)
}

posterior <- function(parm,data){
    return <- list(LP=1,Dev=1,Monitor=1,yhat=model(parm),parm=parm)
}

Sample <- LaplacesDemon(hyperprior, mydata, Initial.Values=c(5), Thinning=10,
                     Iterations=10000, Status=1000,
                     Algorithm="AFSS",
                     Specs=list(A=500, B=NULL, m=100, n=0, w=1)
                     )

plot(Sample,BurnIn=500,mydata,PDF=TRUE,Parms=NULL)

postpredictive <- predict(Sample,posterior,mydata,CPUs=2)

plot(postpredictive, Style="Density",PDF=TRUE)

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


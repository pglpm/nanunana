library('magicaxis')
library('mvtnorm')
library('ellipse')
library('LaplacesDemon')
library('RColorBrewer')

mymu <- 3
mysigma <- 2

draws <- rnorm(10,-3,mysigma)
n <- length(draws)
print(mean(draws))

print(sd(draws)*(n-1)/n)

Data <- list(y=draws, mon.names='predictive distr.', parm.names=c('mu'),N=n)
Data2 <- list(y=100, mon.names='predictive distr.', parm.names=c('mu'),N=n)

Model=function(parm,Data){
    mu.prior <- dnorm(parm, 0, 1000, log=TRUE)
    LL <- sum(dnorm(Data$y,parm,mysigma, log=TRUE))
    LP <- LL + mu.prior
    return=list(LP=LP,Dev=2*LL,Monitor=rnorm(1,parm,mysigma),yhat=rnorm(1,parm,mysigma),parm=parm)}

Fit <- LaplacesDemon(Model, Data, Initial.Values=c(4),
                        Covar=NULL, Iterations=10000, Status=1000, Thinning=1,
                        CPUs=2, Algorithm="AFSS",
                        Specs=list(A=500, B=NULL, m=100, n=0, w=1))


plot(Fit,BurnIn=500,Data,PDF=TRUE,Parms=NULL)

Pred <- predict(Fit,Model,Data2,CPUs=2)

plot(Pred, Style="Density",PDF=TRUE)

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


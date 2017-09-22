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

Data <- list(y=c(10), mon.names='predictive distr.', parm.names=c('mu'),N=1)
Data2 <- list(y=c(20), mon.names='predictive distr.', parm.names=c('mu'),N=1)

Model <- function(parm,Data){
#    mu.prior <- dunif(parm, 0, 1000, log=TRUE)
#    LL <- sum(dnorm(Data$y,parm,mysigma, log=TRUE))
    monitor <- parm
    LP <- dunif(parm,0,1,log=TRUE)
    return=list(LP=1,Dev=1,Monitor=monitor,yhat=Data$y,parm=parm+1)}

Model2 <- function(parm,Data){
#    mu.prior <- dunif(parm, 0, 1000, log=TRUE)
#    LL <- sum(dnorm(Data$y,parm,mysigma, log=TRUE))
    LP <- dunif(parm,0,1,log=TRUE)
    return=list(LP=2,Dev=1,Monitor=1,yhat=parm+Data$y,parm=parm)}

Fit <- LaplacesDemon(Model, Data, Initial.Values=c(4), Thinning=1,
                     Iterations=20, Status=1,
                     Algorithm="AFSS",
                     Specs=list(A=0, B=NULL, m=1, n=0, w=1)
                     )

#plot(Fit,BurnIn=500,Data,PDF=TRUE,Parms=NULL)

Pred11 <- predict(Fit,Model,Data,CPUs=1)
Pred12 <- predict(Fit,Model,Data2,CPUs=1)
Pred21 <- predict(Fit,Model2,Data,CPUs=1)
Pred22 <- predict(Fit,Model2,Data2,CPUs=1)

c(str(Pred11), str(Pred12), str(Pred21), str(Pred22))

stop()

#plot(Pred, Style="Density")

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


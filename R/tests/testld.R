library('magicaxis')
library('mvtnorm')
library('ellipse')
library('LaplacesDemon')

draws <- rnorm(40,1,2)
print(mean(draws))
Data <- list(data=draws, mon.names='', parm.names=c('mu','sigma'),N=length(draws))

Model=function(parm,Data){
    parm[2] <- interval(parm[2],1e-100,100)
    mu.prior <- dnorm(parm[1], 0,1000, log=TRUE)
    sigma.prior <- dhalfcauchy(parm[2], 25, log=TRUE)
    LL <- sum(dnorm(Data$data,parm[1],parm[2], log=TRUE))
    LP <- LL + mu.prior + sigma.prior
    return=list(LP=LP,Dev=2*LL,Monitor=1,yhat=1,parm=parm)}

FitCoin2D <- LaplacesDemon(Model, Data, Initial.Values=c(0,1),
                        Covar=NULL, Iterations=10000, Status=1000, Thinning=1,
                        Algorithm="AFSS",
                        Specs=list(A=500, B=NULL, m=100, n=0, w=1))

png('testplot.png')
magcon(FitCoin2D$Posterior2[9000:10000,1],FitCoin2D$Posterior2[9000:10000,2],
       conlevels=c(0.5,0.68,0.95), lty=c(2,1,3))

title(xlab='mu',ylab='sigma')

polygon(ellipse(cov.rob(FitCoin2D$Posterior2)$cov,centre=FitCoin2D$Summary2[,'Mean'],
                level=c(pnorm(1)-pnorm(-1))), col=hsv(v=0,alpha=0.2),border=NA)

points(1,2,col='blue',pch=4)

points(FitCoin2D$Summary2[1,'Median'],FitCoin2D$Summary2[2,'Median'],
col='red',pch=4)
dev.off()

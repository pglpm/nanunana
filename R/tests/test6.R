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

sigma <- 4

mu <- 4
sigma0 <- 1000
mu0 <- 5

draws <- rnorm(3,mu,sigma)
draws0 <- 6
n <- length(draws)
me <- mean(draws)

sigman <- 1/(n/sigma + 1/sigma0)
mun <- sigman*(mu0/sigma0 + n*me/sigma)

sigmat <- sigman + sigma

PGF <- function(data){
    return(rnormv(1,mu0,sigma0))
}


mydata <- list(data=draws, predict=draws0,PGF=PGF, mon.names=c('P(d0)','y'), parm.names=c(expression('mu')),N=n, y=mean(draws))

mydata0 <- list(data=draws, predict=draws0, PGF=PGF, mon.names=c('P(all data)'), parm.names=c(expression('mu')),N=n, y=mean(draws))

model <- function(parm,number=1){
    return <- rexp(number, rate=parm)
}


hyperprior <- function(parm,data){
    parm.prior <- dnormv(parm,mu0,sigma0, log=T)
    LL <- sum(dnormv(data$data,parm,sigma,log=T))
    LP <- LL+ parm.prior
    return <- list(LP=LP, Dev=-2*LL, Monitor=c(dnormv(data$predict,parm,sigma),rnormv(1,parm,sigma))
                  ,yhat=rnormv(1,parm,sigma),parm=parm)
}

hyperprior0 <- function(parm,data){
    parm.prior <- dnormv(parm,mu0,sigma0, log=T)
    LP <- parm.prior
    return <- list(LP=LP, Dev=-2*LP, Monitor=exp(sum(dnormv(data$data,parm,sigma,log=T))) ,yhat=rnormv(1,parm,sigma),parm=parm)
}

posterior <- function(parm,data){
    LP <- 1
    return <- list(LP=LP,Dev=1,Monitor=1,yhat=rnormv(1,parm,sigma),parm=parm)
}


Sample <- LaplacesDemon(hyperprior, mydata, Initial.Values=c(1), Thinning=2,
                     Iterations=100000, Status=10000,
                     Algorithm="AFSS",
                     Specs=list(A=500, B=NULL, m=100, n=0, w=1)
                     )

plot(Sample,BurnIn=500,mydata,PDF=TRUE,Parms=NULL)


#postpredictive <- predict(Sample,posterior,mydata,CPUs=2)

## png('testnormal_normal.png')
## densplot(postpredictive$yhat,
##     adjust=0.1,
##     main='test normal-normal',
##     xlab=expression('mu'))
## xx <- seq(mun-3*sigmat,mun+3*sigmat,0.1)
## yy <- dnormv(xx,mun,sigmat)
## lines(x=xx,y=yy,col='red')
## dev.off()

png('predictive_probability.png')
exactpdata <- dnormv(mydata$predict,mun,sigmat)
pdata <- mean(Sample$Monitor[,1])
densplot(Sample$Monitor[,1],
    adjust=0.002,
    main='predictive probability for data_0 + uncertainty',
    xlab=paste('P(data_0 = ',draws0,' | data_training) = ',pdata,' (red = est., blue = exact)'),
    ylab='p(P)')
abline(v=pdata,col='red')
abline(v=exactpdata,col='blue')
dev.off()

png('predictive_posterior.png')
densplot(Sample$Monitor[,2],
    adjust=0.1,
    main='predictive posterior distribution',
    xlab=expression('data_0 (blue = exact distr)'),
    ylab=expression('P(data_0 | data_training)'))
xx <- seq(mun-3*sigmat,mun+3*sigmat,0.1)
yy <- dnormv(xx,mun,sigmat)
lines(x=xx,y=yy,col='blue')
dev.off()

Sample0 <- LaplacesDemon(hyperprior0, mydata0, Initial.Values=c(1), Thinning=2,
                     Iterations=100000, Status=10000,
                     Algorithm="AFSS",
                     Specs=list(A=500, B=NULL, m=100, n=0, w=1)
                     )


png('model_probability.png')
exactp <- sqrt(sigma)/sqrt((2*pi*sigma)^n*(n*sigma0+sigma))*exp(
-sum(draws^2)/(2*sigma)-mu0^2/(2*sigma0)+(sigma0*n^2*me^2/sigma + sigma*mu0^2/sigma0 + 2*n*me*mu0)/(2*(n*sigma0+sigma)))
pmodel <- mean(Sample0$Monitor[,1])
densplot(Sample0$Monitor[,1],
    adjust=max(Sample0$Monitor[,1])/50,
    main='probability of model + uncertainty',
    xlab=paste('P(model | data) = ',pmodel,' (red=est., blue=exact)'))
abline(v=pmodel,col='red')
abline(v=exactp,col='blue')
dev.off()

plot(Sample0,BurnIn=500,mydata0,PDF=TRUE,Parms=NULL)



stop()

sum(exp(-(mu0-draws)^2/(2*(sigma+sigma0))))/sqrt(2*pi*(sigma+sigma0))


                                       



len <- length(Sample0$Posterior2)
summ <-0
for(i in 1:len){
summ <- summ + prod(len*dnormv(draws,Sample0$Posterior2[i],sigma))
}
summ/len^length(draws)/len


densplot(Sample$Monitor,
    adjust=1,
    main='test',
    xlab='tau')

xx <- seq(-8,8,0.1)
yy <- dt.scaled(xx,df,mean=mu,sd=sd)
lines(x=xx,y=yy,col='red')


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


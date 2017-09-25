## test script for normal model with Jeffreys prior for variances and unif.
## for correlation - t-walk Monte Carlo
library('Matrix')
library('magicaxis')
library('ellipse')
library('MCMCpack')
library('LaplacesDemon')
library('RColorBrewer')
library('mvtnorm')
library('tmvtnorm')
library('magrittr')
library('bayesplot')
library('R2Cuba')
library('ggplot2')
library('lavaan')

## colourblind-friendly palette
mypurpleblue <- '#4477AA'
myblue <- '#66CCEE'
mygreen <- '#228833'
myyellow <- '#CCBB44'
myred <- '#EE6677'
myredpurple <- '#AA3377'
mygrey <- '#BBBBBB'
palette(c(myblue, myred, mygreen, myyellow, myredpurple, mypurpleblue, mygrey, 'black'))
dev.off()

## string for plots
fname <- 'test10'

## density-plot function
densplot <- function (x,adjust=1,...) { density(x,adjust) %>% plot(.,...)}

set.seed(666)

## parameters to generate normal training data
d <- 3 # num graph quantities
## all graph quantities are in [0,1]
lowb <- rep(0,d)
uppb <- rep(1,d)
mud <- c(0.2,0.3,0.2)
sigmad <- matrix(c(1,2,3,0,0.5,1,0,0,0.1), d)
sigmad <- t(sigmad)+sigmad # symmetrize

## training data with means and scatter matrix
nt <- 10 # num. training data
datat <- signif(rtmvnorm(nt,mud,sigmad,lowb,uppb),2)
meant <- colMeans(datat)
covt <- cov(datat)*(nt-1)/nt

## new datum
datanew <- signif(rtmvnorm(1,mud,sigmad,lowb,uppb),2)

## the parameters are:
## means (d) - uniform distribution (broad normal: curvature helps Monte Carlo)
mean0 <- 0.5
sigma0 <- 10
meanv <- 0
sigmav <- log(1e3)
## variances (d) - uniform in logarithm = Jeffreys dv/v
## correlations (d(d-1)/2) - uniform in [-1,1]
nr <- d*(d-1)/2 # num. correlations
np <- 2*d + nr # tot num. params
parm.names <- c( as.parm.names(list(mu=rep(0,d),logvar=rep(0,d))),
    as.parm.names(list(rho=diag(d)))[upper.tri(diag(d))] )

## function to generate initial values
PGF <- function(data){
    mu <- rnorm(3,mean0,sigma0)
    logvar <- rnorm(3,meanv,sigmav)
    rho <- runif(nr,-1,1)
    return(c(mu,logvar,rho))
}

## model data (barely used)
mydata <- list(data=datat, predict=datanew, PGF=PGF,
               mon.names=c('P(d_new)','d_new[1]','d_new[2]','d_new[3]'),
               parm.names=parm.names,
               N=nt, y=meant)


## model
hyperprior0 <- function(parm,data){
    mu <- interval(parm[1:d],-2,3)
    parm[1:d] <- mu
    logvar <- interval(parm[(d+1):(2*d)],-10,10)
    parm[(d+1):(2*d)] <- logvar
    ## ensure that correlations are between -1 and 1
    rho <- interval(parm[(2*d+1):np],-1,1)
    parm[(2*d+1):np] <- rho
    print(paste(signif(c(mu,logvar,rho),2)))
    ## priors 
    mu.prior <- sum(dnorm(mu,mean0,sigma0,log=T))
    logvar.prior <- sum(dnorm(logvar,meanv,sigmav,log=T))
    ## prior for corr. is uniform
    rho.prior <- sum(dnorm(rho,meanv,sigmav*3,log=T)) # dunif(rho,-1,1,log=T)
    ## data likelihood
    cov <- nearPD(getCov(rho, sds=sqrt(exp(logvar)), diagonal=F))$mat
    LL <- sum(dtmvnorm(data$data, mu, cov, lowb,uppb, log=T))
    LP <- LL + mu.prior + logvar.prior + rho.prior
    return <- list(LP=LP, Dev=-2*LL,
                   ## sample also: p(d=datanew|mu,sigma), and p(d|mu,sigma)
                   Monitor=c(dtmvnorm(data$predict,mu,cov,lowb,uppb),rtmvnorm(1,mu,cov,lowb,uppb,algorithm='gibbs')),
                  yhat=1,parm=parm)
}


## Monte Carlo sampling
initval <- c(rep(mean0,d),rep(meanv,d),rep(meanv,nr))
Sampleinitial <- LaplacesDemon(hyperprior0, mydata, initval,
                        Covar=NULL,
                        Thinning=2,
                        Iterations=50, Status=10,
##                        Algorithm="NUTS", Specs=list(A=1000, delta=0.6, epsilon=1, Lmax=5)
##                        Algorithm="AFSS", Specs=list(A=1000, B=NULL, m=100, n=0, w=1)
                        Algorithm="AFSS", Specs=list(A=100, B=NULL, m=100, n=0, w=1)
##                        Algorithm="twalk", Specs=list(SIV=NULL, n1=4, at=6, aw=1.5)
                        )

stop()

Sample0 <- LaplacesDemon(hyperprior0, mydata, as.initial.values(Sampleinitial),
                        Covar=NULL,
                        Iterations=20000, Status=10000, Thinning=4,
                        Algorithm="AFSS", Specs=list(A=0, B=NULL, m=100, n=0, w=1)
##                        Algorithm="NUTS", Specs=list(A=1000, delta=0.6, epsilon=1, Lmax=5)
##                        Algorithm="twalk", Specs=list(SIV=NULL, n1=4, at=6, aw=1.5)
                        )

### plots:

nsamples <- dim(Sample0$Posterior1)[1]
postdist <- rep(0,nsamples)
for(i in 1:nsamples){
    cov <- getCov(Sample0$Posterior1[i,(2*d+1):np], sds=sqrt(exp(Sample0$Posterior1[i,(d+1):(2*d)])), diagonal=F)
    postdist[i] <- rtmvnorm(1,Sample0$Posterior1[i,1:d],cov,lowb,uppb)
}



## predictive distribution as scatter + marginals
png(paste0('predictive_distr_grid_',fname,'.png'))
mcmc_pairs(Sampleinitial$Posterior1)
dev.off()

## predictive distribution as density
png(paste0('predictive_distr_dens_',fname,'.png'))
magcon(Sample0$Monitor[,2],Sample0$Monitor[,3],# xlim=c(-20,20), ylim=c(-40,40),
       conlevels=c(0.05,0.5,0.95), lty=c(2,1,3),
       imcol=brewer.pal(n=9,name='Blues'))
title(xlab=mydata$mon.names[2],ylab=mydata$mon.names[3],main='predictive distribution (Monte Carlo)')
points(Sample0$Summary2[8,'Mean'],Sample0$Summary2[9,'Mean'],
       col='black',pch=4)
for(i in 1:length(datat[,1])){
points(datat[i,1],datat[i,2], col='#BBBBBB',pch=18)
}
dev.off()

## probability for datanew + 'uncertainty'
png(paste0('prob_datanew_',fname,'.png'))
pnew.mean <- mean(Sample0$Monitor[,1])
densplot(Sample0$Monitor[,1],
    adjust=0.0005,
    main='predictive probability for d_new + uncertainty',
    xlab=paste0('P(d_new = (',signif(datanew[1],2),',',signif(datanew[2],2),') | data_training) = ',signif(pnew.mean,2)),
    ylab='p(P)')
abline(v=pnew.mean,col=myred)
dev.off()


## posterior for the parameters
png(paste0('posterior_parameters_',fname,'.png'))
mcmc_pairs(Sample0$Posterior2)
dev.off()

save.image(file=paste0('normalnormaljeffreys_',fname,'.RData'))

stop()


## value of P(d_new | d_training) via quadrature

hyperpriorc <- function(parm,weight){
    var1 <- exp(parm[3])
    var2 <- exp(parm[4])
    cov <- parm[5]*sqrt(var1*var2)
    varm <- matrix(c(var1,cov,cov,var2),2)
    return(exp(sum(dmvnorm(datat,mean=parm[1:2],sigma=varm, log=T))
               ## + dmvnorm(parm[1:2],mean=mean0,sigma=sigma0, log=T)+
               ## dnorm(parm[3],0,log(1e6),log=T)+
               ## dnorm(parm[4],0,log(1e6),log=T)
               ))
}

normf <- vegas(5,1,hyperpriorc,
      lower=c(-1e6,-1e6,-1e4,-1e4,-1), upper=c(1e6,1e6,1e4,1e4,1),
      rel.tol= 1e-2, abs.tol= 1e-40, max.eval=1000000,
      flags= list(verbose=1, final=0))

hyperpriorcd <- function(parm,weight){
    var1 <- exp(parm[3])
    var2 <- exp(parm[4])
    cov <- parm[5]*sqrt(var1*var2)
    varm <- matrix(c(var1,cov,cov,var2),2)
    return(exp(dmvnorm(datanew,mean=parm[1:2],sigma=varm,log=T)
               + sum(dmvnorm(datat,mean=parm[1:2],sigma=varm, log=T))
               ## + dmvnorm(parm[1:2],mean=mean0,sigma=sigma0, log=T)+
               ## dnorm(parm[3],0,log(1e6),log=T)+
               ## dnorm(parm[4],0,log(1e6),log=T)
               ))
}

normfd <- cuhre(5,1,hyperpriorcd,
      lower=c(-1e6,-1e6,-1e6,-1e6,-1), upper=c(1e6,1e6,1e6,1e6,1),
      rel.tol= 1e-2, abs.tol= 1e-40, max.eval=1000000,
      flags= list(verbose=1, final=0))

pnew.quad <- normfd$value/normf$value



##### garbage and test scripts #####


pnew.exact <- mvtnorm::dmvt(datanew, delta=meant0, sigma=covt0*(nt0+1)/(nt0-d+1), df=nt0-d+1, type='shifted',log=F)

hyperprior <- function(parm,data){
    mu <- parm[1:2]
    parm[3] <- interval(parm[3],0,Inf)
    parm[5] <- interval(parm[5],0,Inf)
    parm[4] <- interval(parm[4]/sqrt(parm[3]*parm[5]),-1,1)*sqrt(parm[3]*parm[5])
    sigma <- matrix(c(parm[c(3,4,4,5)]), 2)
    sigma.prior <- dinvwishart(sigma,n0,scatterm0*n0, log=T)
    mu.prior <- dmvnorm(mu,mean=mean0,sigma=sigma/n0, log=T)
    LL <- sum(dmvnorm(data$data,mean=mu,sigma=sigma,log=T))
    LP <- LL + mu.prior + sigma.prior
    return <- list(LP=LP, Dev=-2*LL,
                   Monitor=c(dmvnorm(data$predict,mean=mu,sigma=sigma),rmvnorm(1,mean=mu,sigma=sigma)),
                  yhat=rmvnorm(1,mean=mu,sigma=sigma),parm=parm)
}

Sample <- LaplacesDemon(hyperprior, mydata, Initial.Values,
                        Covar=NULL,
                        Thinning=10,
                        Iterations=10000, Status=1000,
                        Algorithm="AFSS",
                        Specs=list(A=1000, B=NULL, m=100, n=0, w=1)
                        )


## plot(Sample,BurnIn=500,mydata,PDF=TRUE,Parms=NULL)

stop()
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
pdata <- signif(mean(psamp),2)
densplot(psamp,
    adjust=0.0002,
    main='predictive probability for data_0 + uncertainty',
    xlab=paste0('P(data_0 = (',signif(datanew[1],2),',',signif(datanew[2],2),') | data_training) = ',pdata,' (red = est., blue = exact)'),
    ylab='p(P)')
abline(v=pdata,col='red')
abline(v=pnew.exact0,col='blue', lty=2)
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


#### garbage & temp scripts ####


tsample <- Sample0
for(i in 1:4){
    for(j in (i+1):5){
png(paste0('testplottrajectory',i,j,'.png'))
magplot(tsample$Posterior1[seq(1,dim(tsample$Posterior1)[1],length.out=200),c(i,j)],type='l',col=hsv(alpha=0.3),xlab=mydata$parm.names[i],ylab=mydata$parm.names[j])
dev.off()
    }}

for(i in 1:4){
    for(j in (i+1):5){
png(paste0('testplotdens',i,j,'.png'))
magcon(tsample$Posterior2[,i],tsample$Posterior2[,j],
       conlevels=c(0.5,0.68,0.95), lty=c(2,1,3), imcol=brewer.pal(n=9,name='Blues'))
title(xlab=mydata$parm.names[i],ylab=mydata$parm.names[j])
##polygon(ellipse(cov.rob(tsample$Posterior2)$cov,centre=tsample$Summary2[,'Mean'],
##               level=c(pnorm(1)-pnorm(-1))), col=hsv(v=0,alpha=0.2),border=NA)
##points(mymu,mysigma,col='blue',pch=4)
points(tsample$Summary2[i,'Median'],tsample$Summary2[j,'Median'],
       col='black',pch=4)
dev.off()
    }}

tsample <- Sample0
png(paste0('testplotdensposterior.png'))
magcon(tsample$Monitor[,2],tsample$Monitor[,3], xlim=c(-20,20), ylim=c(-40,40),
       conlevels=c(0.5,0.68,0.95), lty=c(2,1,3), imcol=brewer.pal(n=9,name='Blues'))
title(xlab=mydata$mon.names[2],ylab=mydata$mon.names[3],main='Monte Carlo')
##polygon(ellipse(cov.rob(tsample$Posterior2)$cov,centre=tsample$Summary2[,'Mean'],
##               level=c(pnorm(1)-pnorm(-1))), col=hsv(v=0,alpha=0.2),border=NA)
##points(mymu,mysigma,col='blue',pch=4)
points(tsample$Summary2[8,'Mean'],tsample$Summary2[9,'Mean'],
       col='black',pch=4)
dev.off()






xgrid <- seq(meant0[1]-1*covt0[1,1],meant0[1]+1*covt0[1,1],length.out=20)
ygrid <- seq(meant0[2]-1*covt0[2,2],meant0[2]+1*covt0[2,2],length.out=20)
tgrid <- expand.grid(x=xgrid,y=ygrid)
values <- mvtnorm::dmvt(as.matrix(tgrid),delta=(meant0),sigma=covt0*(nt0+1)/(nt0-d+1), df=nt0-d+1, type='shifted',log=F)
png('contour.png')
contour(x = xgrid, y = ygrid, z = matrix(values, nrow=20, ncol=20),nlevels=50)
dev.off()

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

testPGF <- function(data){ return(interval(rnorm(1,2,1),0,1)) }

testdata <- list(data=0, PGF=testPGF, mon.names=c(' '), parm.names=c('theta'), N=1)

testdistr <- function(parm,data){
   parm <- interval(parm,0,1)
#    LP <- if(parm >0 & parm <1){(dnorm(parm,2,1,log=T))} else {1e-6*dnorm(parm,0.5,10000,log=T)}
    LP <- dnorm(parm,2,1,log=T)
    return <- list(LP=LP, Dev=-2*LP, Monitor=c(1), yhat=1, parm=parm)
}


Initial.Values <- c(0.5)
time2 <- system.time(
    testsample <- LaplacesDemon(testdistr, testdata, Initial.Values,
                                Covar=NULL,
                                Thinning=4,
                                Iterations=10000, Status=5000,
                                Algorithm="AFSS", Specs=list(A=0, B=NULL, m=100, n=0, w=1)
                                )
)

## probability for datanew + 'uncertainty'
testexact <- rtnorm(100000,2,1,0,1)
xx=seq(0,1, length.out=101)
testexact2 <- dtnorm(xx,2,1,0,1)
dat <- data.frame(dens=c(testsample$Posterior2,testsample2$Posterior2,testexact), lines=c(rep('MC',length(testsample$Posterior2)),rep('MC2',length(testsample2$Posterior2)), rep('exact',length(testexact))))
dat2 <- data.frame(yv=testexact2, xv=xx)
png('testdistr.png')
ggplot() + geom_line(data=dat2, aes(x=xv,y=yv)) + geom_density(data=dat, aes(x = dens, fill = lines, alpha =0.5), adjust=1)
dev.off()


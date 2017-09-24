## test script for normal model with normal-Wishart prior (cf minka1998_r2001)
library('pacman')
p_load('magicaxis')
p_load('ellipse')
p_load('MCMCpack')
p_load('LaplacesDemon')
p_load('RColorBrewer')
p_load('mvtnorm')
p_load('magrittr')
p_load('bayesplot')

## density-plot function
densplot <- function (x,adjust=1,...) { density(x,adjust) %>% plot(.,...)}

set.seed(666)

## parameters to generate normal training data
d <- 2
mud <- c(1,2)
sigmad <- matrix(c(3,2,2,5), 2)

## training data with means and scatter matrix
nt <- 10
datat <- signif(rmvnorm(nt,mud,sigmad),2)
meant <- colMeans(datat)
hmeant <- sweep(datat,2,meant,'-')
scattermt <- t(hmeant) %*% hmeant

## reference-prior parameters (mean and variance matrix)
n0 <- 4
mean0 <- c(0,0)
scatterm0 <- diag(c(1,1))

## combination of training data and reference parameters

nt0 <- nt + n0
meant0 <- (nt*meant + n0*mean0)/nt0
scattermt0 <- (nt*scattermt + n0*scatterm0)/nt0 +
    nt*n0*((meant-mean0) %*% t(meant-mean0))/nt0^2

## extra datum
datanew <- signif(rmvnorm(1,mud,sigmad),2)

## exact results:

## P(datanew | datat, prior) is t-student
pnew.exact0 <- mvtnorm::dmvt(datanew, delta=meant, sigma=scattermt*(nt+1)/(nt-d+1), df=nt-d+1, type='shifted',log=F)

## function to generate initial values (don't know how it works)
PGF <- function(data){
    sigma <- rinvwishart(n0,scatterm0*n0)
    mu <- rmvnorm(1,mean=mean0,sigma=sigma/n0)
    return(c(mu,c(sigma)[c(1,2,4)]))
}

## model data - barely used
mydata <- list(data=datat, predict=datanew, PGF=PGF,
               mon.names=c('P(d_new)','d_new[1]','d_new[2]'),
               parm.names=c('mu[1]','mu[2]','sigma[1,1]','sigma[1,2]','sigma[2,2]'),
               N=nt, y=meant)


## model
hyperprior0 <- function(parm,data){
    mu <- parm[1:2]
    ## ensure that variance matrix is positive-definite
    parm[3] <- interval(parm[3],1e-6,1e6)
    parm[5] <- interval(parm[5],1e-6,1e6)
    parm[4] <- interval(parm[4]/sqrt(parm[3]*parm[5]),-1,1)*sqrt(parm[3]*parm[5])
    sigma <- matrix(c(parm[c(3,4,4,5)]), 2)
    ## prior for variance is inverse Wishart
    sigma.prior <- dinvwishart(sigma,nt,scattermt*nt, log=T)
    ## prior for mean is normal
    mu.prior <- dmvnorm(mu,mean=meant,sigma=sigma/nt, log=T)
    LP <- mu.prior + sigma.prior
    return <- list(LP=LP, Dev=-2*LP,
                   ## sample also: p(d=datanew|mu,sigma), and p(d|mu,sigma)
                   Monitor=c(dmvnorm(datanew,mean=mu,sigma=sigma),rmvnorm(1,mean=mu,sigma=sigma)),
                  yhat=rmvnorm(1,mean=mu,sigma=sigma),parm=parm)
}

## Monte Carlo sampling
Initial.Values <- c(0,0,1,0,1)
Sample0 <- LaplacesDemon(hyperprior0, mydata, Initial.Values,
                        Covar=NULL,
                        Thinning=10,
                        Iterations=10000, Status=1000,
##                        Algorithm="NUTS", Specs=list(A=1000, delta=0.6, epsilon=1, Lmax=5)
                        Algorithm="AFSS", Specs=list(A=1000, B=NULL, m=100, n=0, w=1)
                        )

### plots:

## predictive distribution as scatter + marginals
png('predictive_distr_grid.png')
mcmc_pairs(Sample0$Monitor[,2:3])
dev.off()

## predictive distribution as density
png('predictive_distr_dens.png')
magcon(Sample0$Monitor[,2],Sample0$Monitor[,3], xlim=c(-20,20), ylim=c(-40,40),
       conlevels=c(0.05,0.5,0.95), lty=c(2,1,3),
       imcol=brewer.pal(n=9,name='Blues'))
title(xlab=mydata$mon.names[2],ylab=mydata$mon.names[3],main='predictive distribution (Monte Carlo)')
points(Sample0$Summary2[8,'Mean'],Sample0$Summary2[9,'Mean'],
       col='black',pch=4)
for(i in 1:length(datat[,1])){
points(datat[i,1],datat[i,2], col='gray',pch=5)
}
dev.off()

## probability for datanew + 'uncertainty'
png('prob_datanew.png')
pnew.mean <- mean(Sample0$Monitor[,1])
densplot(Sample0$Monitor[,1],
    adjust=0.0001,
    main='predictive probability for d_new + uncertainty',
    xlab=paste0('P(d_new = (',signif(datanew[1],2),',',signif(datanew[2],2),') | data_training) = ',signif(pnew.mean,2),' (red = est., blue = exact)'),
    ylab='p(P)')
abline(v=pnew.mean,col='red')
abline(v=pnew.exact0,col='blue', lty=2)
dev.off()


## posterior for the parameters
png('posterior_parameters.png')
mcmc_pairs(Sample0$Posterior2,pch=4)
dev.off()


## exact predictive distribution as density
sampleexact <- mvtnorm::rmvt(100000, delta=meant, sigma=scattermt*(nt+1)/(nt-d+1), df=nt-d+1, type='shifted')
png(paste0('predictive_distr_dens_exact.png'))
magcon(sampleexact[,1],sampleexact[,2], xlim=c(-20,20), ylim=c(-40,40),
       conlevels=c(0.05,0.5,0.95), lty=c(2,1,3),
       imcol=brewer.pal(n=9,name='Blues'))
title(xlab=mydata$mon.names[2],ylab=mydata$mon.names[3], main='predictive distribution (exact)')
points(meant[1],meant[2], col='black',pch=4)
for(i in 1:length(datat[,1])){
points(datat[i,1],datat[i,2], col='gray',pch=5)
}
dev.off()


stop()





##### garbage and test scripts #####


pnew.exact <- mvtnorm::dmvt(datanew, delta=meant0, sigma=scattermt0*(nt0+1)/(nt0-d+1), df=nt0-d+1, type='shifted',log=F)

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






xgrid <- seq(meant0[1]-1*scattermt0[1,1],meant0[1]+1*scattermt0[1,1],length.out=20)
ygrid <- seq(meant0[2]-1*scattermt0[2,2],meant0[2]+1*scattermt0[2,2],length.out=20)
tgrid <- expand.grid(x=xgrid,y=ygrid)
values <- mvtnorm::dmvt(as.matrix(tgrid),delta=(meant0),sigma=scattermt0*(nt0+1)/(nt0-d+1), df=nt0-d+1, type='shifted',log=F)
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


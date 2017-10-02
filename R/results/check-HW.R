## test script for normal model 
library('pacman')
library('magicaxis')
library('ellipse')
library('MCMCpack')
library('LaplacesDemon')
library('RColorBrewer')
library('mvtnorm')
library('magrittr')
library('bayesplot')
library('ggplot2')
mypurpleblue <- '#4477AA'
myblue <- '#66CCEE'
mygreen <- '#228833'
myyellow <- '#CCBB44'
myred <- '#EE6677'
myredpurple <- '#AA3377'
mygrey <- '#BBBBBB'
palette(c(myblue, myred, mygreen, myyellow, myredpurple, mypurpleblue, mygrey, 'black'))
dev.off()

## base filename to save results and title for plots
filename <- 'logit-HW_healthy'
ptitle <- 'healthy (logit-normal + Huang-Wang)'
datafile <- 'data_H.dat'

## seed for random generator
set.seed(666)


## data & data names
y.names <- scan('graph_quantities.dat', what="character", sep=",")[1:3]
datam <- t(read.matrix(datafile))[,1:3]
colnames(datam) <- y.names
rownames(datam) <- sprintf("id[%d]",seq(1:dim(datam)[1]))
## logit of data
data <- logit(datam)
d <- dim(data)[2] # num parms
N <- dim(data)[1] # num inds
dpos <- choose(1:d +1, 2) # (1:d)*((1:d)+1)/2  position of diagonal elements
nr <- d*(d-1)/2 # num. correlations
np <- 2*d + nr # total num. parameters
d1 <- d+1
d2 <- 2*d
d3 <- 2*d+1
## suff statistics: logit-means, -stds, -corrs
dmean <- colMeans(data)
dcov <- cov(data)*(N-1)/N
dstd <- sqrt(diag(dcov))
dcor <- Cov2Cor(dcov)[upper.tri(dcov,diag=F)]

## new datum - we take the mean of the data we have just for check
ynew <- dmean

## hyperparameters for hyperprior
meanmu <- 0 # mean for mu
sigmamu <- 10^2 # variance for mu
shiftsd <- 6 # shift for sd
Ahw <- rep(1e6,d) # scale hyperparam for HuangWand distribution
ahw <- rep(1,d) # scale param for HuangWand distribution
nuhw <- 2 # d.f. for HuangWand; 2 ensures uniform marginals for corrs

## parameters:  means + elements of Cholesky decomposition
parm.names <- as.parm.names(list(mu=rep(0,d), lnU=rep(0,d), U=matrix(0,d,d),gamma=rep(0,d)), uppertri=c(0,0,1,0))
parm.names <- parm.names[-(d2+dpos)]

# parameters we want to see
tparm.names <- as.parm.names(list(mu=rep(0,d), sigma=rep(0,d), rho=matrix(0,d,d)), uppertri=c(0,0,1))
tparm.names <- tparm.names[-(d2+dpos)]

## quantities to be monitored during Monte Carlo: P(d_new), posterior distr
mon.names <- c(y.names,'p(d_new)')

## function to generate initial values (don't know how it works)
PGF <- function(data){
    mu <- rnorm(d,meanmu,sigmamu)
    U <- rhuangwandc(nu=nuhw,a=ahw,A=Ahw)
    gamma <- runif(d)
    return(c(mu,log(diag(U)),U[upper.tri(U)],gamma))
}
Initial.Values <- c(rep(meanmu,d),rep(0,d),rep(0,nr))

## model data, input to the Monte Carlo algorithm
mydata <- list(y=data, PGF=PGF,
               parm.names=parm.names,
               mon.names=mon.names,
#               pos.mu=pos.mu, pos.sigma=pos.sigma, pos.rho=pos.rho,
               N=N, yhat=ynew)

## model
hyperprior0 <- function(parm,data){
    ## prior for mean is normal
    mu <- parm[1:d]
    mu.prior <- sum(dnorm(mu,meanmu,sigmamu, log=T))
    ## prior for Cholesky matrix
    U <- diag(exp(parm[d1:d2]))
    U[upper.tri(U)] <- parm[d3:np]
    gamma <- interval(parm[(np+1):(np+d)], 1e-100,Inf)
    parm[(np+1):(np+d)] <- gamma
    U.prior <- dhuangwandc(U,nu=nuhw,a=gamma,A=Ahw,log=T) 
    ## construct covariance matrix for data log-likelihood
    #print(U) # debug
    covm <- t(U) %*% U
    ## log-likelihood
    LL <- 0 # sum(dmvnorm(data$y, mu, covm, log=T))
    ## log-posterior
    LP <- LL + mu.prior + U.prior
    return <- list(LP=LP, Dev=-2*LL,
                   ## monitor posterior and probability uncertainty for new datum
                   Monitor=c(invlogit(rmvnorm(1,mu,covm)),
                             dmvnorm(data$yhat,mu,covm)),
                   yhat=1,
                   parm=parm)
}


## Monte Carlo sampling:

## First short adaptive sampling
sampleinitial <- LaplacesDemon(hyperprior0, mydata, Initial.Values,
                        Covar=NULL,
                        Thinning=1,
                        Iterations=1000, Status=100,
##                        Algorithm="NUTS", Specs=list(A=1000, delta=0.6, epsilon=1, Lmax=5)
                        Algorithm="AFSS", Specs=list(A=500, B=NULL, m=100, n=0, w=1)
                        )

## Longer sampling, with two chains on parallel CPUs
sample2 <- LaplacesDemon.hpc(hyperprior0, mydata, as.initial.values(sampleinitial),
                        Covar=sampleinitial$Covar,
                        Thinning=1,
                        Iterations=5e5, Status=1000,
                        Chains=2,CPUs=2,Packages=c('mvtnorm'),LogFile=paste0(filename,'_LD_log'),#Type="MPI",
##                        Algorithm="NUTS", Specs=list(A=1000, delta=0.6, epsilon=1, Lmax=5)
                        Algorithm="AFSS", Specs=list(A=0, B=NULL, m=100, n=0, w=1)
                        )

## combine the chains
sample0 <- Combine(sample2,mydata)


## transform hyperparameters in understandable ones: stds, corrs
samples <- sample0$Posterior1
nsamples <- dim(samples)[1]
for(i in 1:nsamples){
    covm <- diag(exp(samples[i,d1:d2]))
    covm[upper.tri(covm)] <- samples[i,d3:np]
    covm <- t(covm) %*% covm
    samples[i,d1:d2] <- sqrt(diag(covm))
    samples[i,d3:np] <- Cov2Cor(covm)[upper.tri(covm)]
}

## Save all data and results
save.image(file=paste0(filename,'.RData'))

### plots:

## posterior predictive for the quantities
binningdiv <- 10
pdf(paste0('predictive_posterior_',filename,'.pdf'))
for(i in 1:d){
    posterior <- sample0$Monitor[,i]
    dat <- data.frame(x=posterior)
    print(
        ggplot() + geom_histogram(data=dat,aes(x=x,y=..density..),
                                  binwidth=1/50) +
        labs(x=mon.names[i],title=ptitle) +
        geom_vline(xintercept=mean(posterior),colour=myred)
        + xlim(0,1)
        #+geom_vline(xintercept=median(posterior),colour=mygreen)
    )
}
for(j in 1:(d-1)){
    for(i in (j+1):d){
        posterior <- sample0$Monitor[,j]
        posterior2 <- sample0$Monitor[,i]
        dat <- data.frame(x=posterior,y=posterior2)
        datpoints <- data.frame(x=datam[,j],y=datam[,i])
        print(
            ggplot() + stat_bin2d(data=dat,aes(x=x,y=y,fill=..density..),
                                  #binwidth=rep(1/25,2),drop=F
                                  breaks=seq(0,1,length.out=25) #, trans="log10"
                                  ) +
            scale_fill_gradientn(colours=brewer.pal(n=9,name='Blues')) +
            labs(x=mon.names[j],y=mon.names[i],title=ptitle) +
             xlim(0,1) + ylim(0,1) +
            geom_point(data=data.frame(x=mean(posterior),y=mean(posterior2))
                      ,aes(x,y),colour=myred,size=4,shape=16)
            + geom_point(data=datpoints,aes(x=x,y=y),
                         colour=myyellow,size=4,shape=5)
        )
    }}
dev.off()

## probability for new datum + its "uncertainty"
posterior <- sample0$Monitor[-(1:d)]
pdf(paste0('prob_newdatum_',filename,'.pdf'))
    dat <- data.frame(x=posterior)
    print(
        ggplot() + geom_histogram(data=dat,aes(x=x,y=..density..),
                                  binwidth=sd(posterior)/binningdiv) +
        labs(x=paste0('p(d_new) = ',signif(mean(posterior),2)),title=ptitle) +
        geom_vline(xintercept=mean(posterior),colour=myred)
        #+geom_vline(xintercept=median(posterior),colour=mygreen)
    )
dev.off()

## posterior for the parameters
binningdiv <- 10
pdf(paste0('hyperparameter_posterior_',filename,'.pdf'))
for(i in 1:np){
    posterior <- samples[,i]
    dat <- data.frame(x=posterior)
    print(
        ggplot() + geom_histogram(data=dat,aes(x=x,y=..density..),
                                  binwidth=sd(posterior)/binningdiv) +
        labs(x=tparm.names[i],title=ptitle) +
        geom_vline(xintercept=mean(posterior),colour=myred)
        #+geom_vline(xintercept=median(posterior),colour=mygreen)
    )
}
for(j in 1:(np-1)){
    for(i in (j+1):np){
        ## if there are more than 10 hyperparameters we choose a random subset of marginals to show
        if(np<11 || (np>10 && runif(1)<50/choose(np,2))){
        posterior <- samples[,j]
        posterior2 <- samples[,i]
        dat <- data.frame(x=posterior,y=posterior2)
        print(
            ggplot() + stat_bin2d(data=dat,aes(x=x,y=y,fill=..density..),
                                  binwidth=2*c(sd(posterior),
                                               sd(posterior2))/binningdiv
                                 #, trans="log10"
                                  ) +
            scale_fill_gradientn(colours=brewer.pal(n=9,name='Blues')) +
            labs(x=tparm.names[j],y=tparm.names[i],title=ptitle) +
            geom_point(data=data.frame(x=mean(posterior),y=mean(posterior2))
                      ,aes(x,y),colour=myred,size=4,shape=16)
        )
        }
    }}
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

sample0 <- LaplacesDemon(hyperprior0, mydata0, Initial.Values=c(1), Thinning=2,
                     Iterations=100000, Status=10000,
                     Algorithm="AFSS",
                     Specs=list(A=500, B=NULL, m=100, n=0, w=1)
                     )


png('model_probability.png')
exactp <- sqrt(sigma)/sqrt((2*pi*sigma)^n*(n*sigma0+sigma))*exp(
-sum(draws^2)/(2*sigma)-mu0^2/(2*sigma0)+(sigma0*n^2*me^2/sigma + sigma*mu0^2/sigma0 + 2*n*me*mu0)/(2*(n*sigma0+sigma)))
pmodel <- mean(sample0$Monitor[,1])
densplot(sample0$Monitor[,1],
    adjust=max(sample0$Monitor[,1])/50,
    main='probability of model + uncertainty',
    xlab=paste('P(model | data) = ',pmodel,' (red=est., blue=exact)'))
abline(v=pmodel,col='red')
abline(v=exactp,col='blue')
dev.off()

plot(sample0,BurnIn=500,mydata0,PDF=TRUE,Parms=NULL)



stop()


#### garbage & temp scripts ####


tsample <- sample0
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

tsample <- sample0
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


                                       



len <- length(sample0$Posterior2)
summ <-0
for(i in 1:len){
summ <- summ + prod(len*dnormv(draws,sample0$Posterior2[i],sigma))
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


system.time(
    for(i in 1:5e6){
        cho <- totl
        cho[dpos] <- exp(cho[dpos])
        matr <- matrix(0,d,d)
        matr[upper.tri(matr,diag=T)] <- cho
    }
)

totl <- 1:10
totd <- totl[dpos]
toto <- totl[-dpos]

system.time(
    for(i in 1:5e6){
        matr <- diag(exp(totl[dpos]))
        matr[upper.tri(matr)] <- totl[-dpos]
    }
) #faster





crv <- rlogis(6)
rv <- 2*invlogit(crv)-1


prior.lrho <- sum(dlogis(crv,log=T))
rho <- diag(d)
rho[lower.tri(rho)] <- rv

rho <- as.symmetric.matrix(rho)

system.time(for(i in 1:1e5){
                test <- 1:10
                rho <- as.symmetric.matrix(1:10)
            })

system.time(for(i in 1:1e5){
rho <- diag(d)
rho[lower.tri(rho)] <- rv
rho[upper.tri(rho)] <- rv
            })


system.time(for(i in 1:1e5){
rho <- rep(1,10)
rho[c(2,3,6,4,7,9)] <- rv
rho <- as.symmetric.matrix(rho)
            })

system.time(for(i in 1:1e5){
rho <- diag(d)
rho[lower.tri(rho)] <- rv
rho[upper.tri(rho)] <- rv
            })

rho[lower.tri(rho),upper.tri(rho)] <- rv

system.time(for(i in 1:1e3){tsum <- sum(dnorm(log(lvtest),0,10,log=T))})

system.time(for(i in 1:1e3){tsum <- sum(dlnorm(lvtest,0,10,log=T))})

dm <- 6
dnp <- dm*(dm+1)/2
nsamp <- 1e5
testq <- matrix(0,nsamp,dnp)
for(i in 1:nsamp){
    matr <- rinvwishart(dm+1,diag(dm))
    testq[i,1:dm] <- log(diag(matr))
    testq[i,(dm+1):dnp] <- Cov2Cor(matr)[upper.tri(matr)]
}
pdf(paste0('test_invwishart_dplus1.pdf'))
    for(i in 1:dnp){
    posterior <- testq[,i]
    densplot(posterior,
    adjust=sd(posterior)/10,
   #  main='predictive probability for d_new + uncertainty',
    xlab=toString(i), ylab='density')
    abline(v=mean(posterior),col=myred)
}
for(j in 1:(dnp-1)){
    for(i in (j+1):dnp){
        posterior <- testq[,j]
        posterior2 <- testq[,i]
        magcon(posterior,posterior2,# xlim=c(-20,20), ylim=c(-40,40),
       conlevels=c(0.05,0.5,0.95), lty=c(2,1,3),
       imcol=brewer.pal(n=9,name='Blues'))
        title(xlab=toString(j),ylab=toString(i))
points(mean(posterior),mean(posterior2),
       col=myred,pch=4)
    }}
dev.off()

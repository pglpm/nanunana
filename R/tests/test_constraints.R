## test script for normal model 
library('pacman')
library('magicaxis')
library('ellipse')
library('MCMCpack')
library('LaplacesDemon')
library('RColorBrewer')
library('mvtnorm')
library('magrittr')
library('ggplot2')
library('bayesplot')

filename <- 'constr_test'
# number of parameters d is below

mypurpleblue <- '#4477AA'
myblue <- '#66CCEE'
mygreen <- '#228833'
myyellow <- '#CCBB44'
myred <- '#EE6677'
myredpurple <- '#AA3377'
mygrey <- '#BBBBBB'
palette(c(myblue, myred, mygreen, myyellow, myredpurple, mypurpleblue, mygrey, 'black'))
dev.off()

## density-plot function
densplot <- function (x,adjust=1,...) { density(x,adjust) %>% plot(.,...)}

# seed for random generator
set.seed(666)



## function to generate initial values (don't know how it works)
PGF <- function(data){
#    mu <- rnorm(d,mean0,sigma0)
#    lvar <- rnorm(d,0,sigmav)
#    lrho <- logit((runif(nr,-1,1)+1)/2)
    return(1)
}


parm.names <- c('1','2')
## model data, input to the Monte Carlo algorithm
mydata <- list(data=1, 
               mon.names=c(''),
               parm.names=parm.names,
               N=1, y=1)


## model
#Initial.Values <- c(rep(0,d),rep(0,nr))
hyperprior0 <- function(parm,data){
    if(parm[1]>0 && parm[2] > 0 && sum(parm)<1){
        LP <- 0} else {
                   LP <- -1e16}
#    parm[1] <- interval(parm[1],0,1)
#    parm[2] <- interval(parm[2],0,1-parm[1])
#    LP <- dmvnorm(parm,c(1/3,1/3),1000*diag(2),log=T)
    return <- list(LP=LP, Dev=-2*LP,
                   ## sample also: p(d=datanew|mu,sigma), and p(d|mu,sigma)
                   Monitor=1,
                  yhat=1,parm=parm)
}

hyperprior2 <- function(parm,data){
    parm[1] <- interval(parm[1],0,1)
    parm[2] <- interval(parm[2],0,1-parm[1])
    if(parm[1]>0 && parm[2] > 0 && sum(parm)<1){
   #     print('** OK **')
        LP <- 0} else {
    #               print('**** BAD ****')
     #              print(parm)
                   LP <- -Inf}

                                        #    LP <- dmvnorm(parm,c(1/3,1/3),1000*diag(2),log=T)
    return <- list(LP=LP, Dev=-2*LP,
                   ## sample also: p(d=datanew|mu,sigma), and p(d|mu,sigma)
                   Monitor=1,
                  yhat=1,parm=parm)
}

hyperprior1 <- function(parm,data){
    if(parm[1]>0 && parm[2] > 0 && sum(parm)<1){
   #     print('** OK **')
        flag <- 0.5} else {
    #               print('**** BAD ****')
     #              print(parm)
                   flag <- -1}
   # parm[1] <- interval(parm[1],0,1)
   # parm[2] <- interval(parm[2],0,1-parm[1])
                                        #    LP <- dmvnorm(parm,c(1/3,1/3),1000*diag(2),log=T)
    LP <- dunif(flag,log=T)
    return <- list(LP=LP, Dev=-2*LP,
                   ## sample also: p(d=datanew|mu,sigma), and p(d|mu,sigma)
                   Monitor=1,
                  yhat=1,parm=parm)
}


## Monte Carlo sampling:


## Longer sampling
sample0 <- LaplacesDemon(hyperprior0, mydata, c(0.1,0.1),
                        Covar=NULL,
                        Thinning=1,
                        Iterations=10000, Status=1000,
##                        Algorithm="NUTS", Specs=list(A=1000, delta=0.6, epsilon=1, Lmax=5)
                        Algorithm="AFSS", Specs=list(A=100, B=NULL, m=100, n=0, w=1)
                        )

png(paste0('posterior_parameterstest_',filename,'.png'))
mcmc_pairs(sample0$Posterior2)
dev.off()


posterior <- sample0$Posterior2
pdf(paste0('posterior_parameters_',filename,'.pdf'))
dm <- length(parm.names)
    for(i in 1:dm){
    posterior <- sample0$Posterior2[,i]
    densplot(posterior,
    adjust=sd(posterior)/10,
   #  main='predictive probability for d_new + uncertainty',
    xlab=parm.names[i], ylab='density')
    abline(v=sample0$Summary2[i,'Mean'],col=myred)
#    abline(v=truevalues[i],col=myblue)
}
for(j in 1:(dm-1)){
    for(i in (j+1):dm){
        posterior <- sample0$Posterior2[,j]
        posterior2 <- sample0$Posterior2[,i]
        dat <- data.frame(x=posterior,y=posterior2)
        print(ggplot() + stat_bin2d(data=dat,aes(x=x,y=y),bins=50) +
            scale_fill_gradientn(colours=brewer.pal(n=9,name='Blues')) +
            labs(x=parm.names[j],y=parm.names[i]))
    }}
dev.off()

dat <- data.frame(x= sample0$Posterior2[,1],y= sample0$Posterior2[,2])
pdf('justatest.pdf')
ggplot() + stat_bin2d(data=dat,aes(x=x,y=y),bins=50) +
    scale_fill_gradientn(colours=brewer.pal(n=9,name='Blues'))
dev.off()

stop()

## Save all data and results
save.image(file=paste0(filename,'.RData'))


### plots:



## posterior for the parameters
if(np<5){
png(paste0('posterior_parameterstest_',filename,'.png'))
mcmc_pairs(sample0$Posterior2)
dev.off()}
##
posterior <- sample0$Posterior2
rlist <- list(rmu=1:d, rvar=(d+1):(2*d), rrho=(2*d+1):np)
pdf(paste0('posterior_parameters_',filename,'.pdf'))
for(j in c('rmu','rvar','rrho')){
    for(i in rlist[[j]]){
    posterior <- sample0$Posterior2[,i]
    densplot(posterior,
    adjust=sd(posterior)/10,
   #  main='predictive probability for d_new + uncertainty',
    xlab=parm.names[i], ylab='density')
    abline(v=sample0$Summary2[i,'Mean'],col=myred)
    abline(v=truevalues[i],col=myblue)
}}
for(j in 1:(np-1)){
    for(i in (j+1):np){
        posterior <- sample0$Posterior2[,j]
        posterior2 <- sample0$Posterior2[,i]
        magcon(posterior,posterior2,# xlim=c(-20,20), ylim=c(-40,40),
       conlevels=c(0.05,0.5,0.95), lty=c(2,1,3),
       imcol=brewer.pal(n=9,name='Blues'))
        title(xlab=parm.names[j],ylab=parm.names[i])
points(sample0$Summary2[j,'Mean'],sample0$Summary2[i,'Mean'],
       col=myred,pch=4)
points(truevalues[j],truevalues[i],
       col=mygreen,pch=4)
    }}
dev.off()

## predictive distribution as scatter + marginals
png(paste0('predictive_distr_grid_',filename,'.png'))
mcmc_pairs(sample0$Monitor[,-1])
dev.off()
##
pdf(paste0('predictive_distr_grid_',filename,'.pdf'))
mcmc_pairs(sample0$Monitor[,-1])
dev.off()

## predictive distribution as density
png(paste0('predictive_distr_dens_',filename,'.png'))
magcon(sample0$Monitor[,2],sample0$Monitor[,3],# xlim=c(-20,20), ylim=c(-40,40),
       conlevels=c(0.05,0.5,0.95), lty=c(2,1,3),
       imcol=brewer.pal(n=9,name='Blues'))
title(xlab=mydata$mon.names[2],ylab=mydata$mon.names[3],main='predictive distribution (Monte Carlo)')
points(sample0$Summary2[8,'Mean'],sample0$Summary2[9,'Mean'],
       col='black',pch=4)
for(i in 1:length(datat[,1])){
points(datat[i,1],datat[i,2], col='#BBBBBB',pch=18)
}
dev.off()
##
pdf(paste0('predictive_distr_dens_',filename,'.pdf'))
magcon(sample0$Monitor[,2],sample0$Monitor[,3],# xlim=c(-20,20), ylim=c(-40,40),
       conlevels=c(0.05,0.5,0.95), lty=c(2,1,3),
       imcol=brewer.pal(n=9,name='Blues'))
title(xlab=mydata$mon.names[2],ylab=mydata$mon.names[3],main='predictive distribution (Monte Carlo)')
points(sample0$Summary2[8,'Mean'],sample0$Summary2[9,'Mean'],
       col='black',pch=4)
for(i in 1:length(datat[,1])){
points(datat[i,1],datat[i,2], col='#BBBBBB',pch=18)
}
dev.off()

## probability for datanew + 'uncertainty'
uncert <- sample0$Monitor[,1]
pnew.mean <- mean(uncert)
png(paste0('prob_datanew_',filename,'.png'))
densplot(uncert,
    adjust=sd(uncert)/10,
    main='predictive probability for d_new + uncertainty',
    xlab=paste0('P(d_new = (',toString(datanew),') | data_training) = ',signif(pnew.mean,2),' (true: ',signif(truepnew,2),')'),
    ylab='p(P)')
abline(v=pnew.mean,col=myred)
abline(v=truepnew,col=myblue)
dev.off()
##
pdf(paste0('prob_datanew_',filename,'.pdf'))
densplot(uncert,
    adjust=sd(uncert)/10,
    main='predictive probability for d_new + uncertainty',
    xlab=paste0('P(d_new = (',toString(datanew),') | data_training) = ',signif(pnew.mean,2),' (true: ',signif(truepnew,2),')'),
    ylab='p(P)')
abline(v=pnew.mean,col=myred)
abline(v=truepnew,col=myblue)
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

parmn <- as.parm.names(list(logvariance=rep(0,dm), rho=diag(dm)),uppertri=c(0,1))
parmn <- parmn[-(dm+choose(1:dm +1, 2))]

dm <- 6
dnp <- dm*(dm+1)/2
nsamp <- 1e5
testq <- matrix(0,nsamp,dnp)
for(i in 1:nsamp){
    matr <- rhuangwand(nu=2,a=rep(1,dm),A=rep(100,dm))
    testq[i,1:dm] <- sqrt(diag(matr))-6
    testq[i,(dm+1):dnp] <- Cov2Cor(matr)[upper.tri(matr)]
}
pdf('test_loghuangwand.pdf')
    for(i in 1:dm){
        posterior <- testq[,i]
    densplot(posterior,
    adjust=0.2,xlim=c(-6,6),
   #  main='predictive probability for d_new + uncertainty',
    xlab=parmn[i], ylab='p')
    abline(v=mean(posterior),col=myred)
    }
    for(i in (dm+1):dnp){
        posterior <- testq[,i]
    densplot(posterior,
    adjust=sd(posterior)/10,
   #  main='predictive probability for d_new + uncertainty',
    xlab=parmn[i], ylab='p')
    abline(v=mean(posterior),col=myred)
}
for(j in 1:(dnp-1)){
    for(i in (j+1):dnp){
        posterior <- testq[,j]
        posterior2 <- testq[,i]
        magcon(posterior,posterior2,# xlim=c(-20,20), ylim=c(-40,40),
       conlevels=c(0.05,0.5,0.95), lty=c(2,1,3),
       imcol=brewer.pal(n=9,name='Blues'))
        title(xlab=parmn[j],ylab=parmn[i])
points(mean(posterior),mean(posterior2),
       col=myred,pch=4)
    }}
dev.off()

## test script for normal model 
library('MCMCpack')
library('LaplacesDemon')
library('RColorBrewer')
library('bayesplot')
library('ggplot2')
library('Matrix')
library('Gmedian')
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
filename <- 'logit-unif_schizo_all_rdmh'
ptitle <- 'logit-uniform model, Schizo all (RDMH)'
datafile <- 'data_S.dat'

## seed for random generator
set.seed(999)

## data & data names
d <- dim(read.matrix(datafile))[1] # num. graph quant.
N <- dim(read.matrix(datafile))[2] # num. indiv.
y.names <- scan('graph_quantities.dat', what="character", sep=",")[1:d]
datam <- matrix(0, N, d)
datam[,] <- t(read.matrix(datafile))[1:N,1:d]
colnames(datam) <- y.names
rownames(datam) <- sprintf("id[%d]",1:N)
## logit of data
data <- logit(datam)
d1 <- d+1
d2 <- 2*d
d3 <- 2*d+1
nr <- d*(d-1)/2 # num. correlations
np <- 2*d + nr # total num. parameters
dnp <- 3*d + nr # total num. parameters
dpos <- choose((1:d) +1, 2) # (1:d)*((1:d)+1)/2  position of diagonal elements
odpos <- (1:(d+nr))[-dpos] # off-diag
ddpos <- d2+dpos
## suff statistics: logit-means, -stds, -corrs
dmean <- colMeans(data)
dcov <- cov(data)*(N-1)/N
dstd <- sqrt(diag(dcov))
#dcor <- Cov2Cor(dcov)[upper.tri(dcov,diag=F)]

## new datum - we take the mean of the data we have just for check
ynew <- dmean

## hyperparameters for hyperprior
meanmu <- 0 # mean for mu
stdmu <- 10 # std for mu
lsigmaa <- -5 # lower bound logsigma
lsigmab <- 5 # upper bound logsigma
nu <- d+1 # scale hyperparam for inv-Wishart
S <- diag(d) # scale matrix for inv-Wishart

## parameters:  means + elements of Cholesky decomposition
pmu <- 1:d
pls <- d1:d2
pU <- d3:dnp
pdU <- d2+dpos
podU <- d2+odpos
parm.names <- as.parm.names(list(mu=rep(0,d), sigma=rep(0,d), U=matrix(0,d,d)), uppertri=c(0,0,1))
temp <- matrix('',d,d)
temp[upper.tri(temp,diag=T)] <- parm.names[pU]
parm.names[pU] <- lower.triangle(t(temp),diag=T)
##
temp <- matrix(0,d,d)
temp[lower.tri(temp,diag=T)] <- pU
pdU <- diag(temp)
podU <-  lower.triangle(temp)
#parm.names <- parm.names[-(d2+dpos)

# parameters we want to see
tparm.names <- as.parm.names(list(mu=rep(0,d), sigma=rep(0,d), rho=matrix(0,d,d)), uppertri=c(0,0,1))
tparm.names <- c(y.names, tparm.names[-(d2+dpos)])



mon.names <- c('')

## function to generate initial values
## PGF2 <- function(data){
##     mu <- rnorm(d,meanmu,stdmu)
##     lsigma <- runif(d,lsigmaa,lsigmab)
##     U <- rinvwishart(d+1,diag(d))
##     return(c(mu, lsigma, lower.triangle(U,diag=T)))
## }
PGF <- function(data){
    mu <- rnorm(d,dmean,sd(dmean))
    lsigma <- interval(rnorm(d,log(dstd),sd(log(dstd))),lsigmaa,lsigmab)
    U <- rinvwishart(d+1,diag(d))
    return(c(mu, lsigma, lower.triangle(U,diag=T)))
}

mydata <- list(y=data, PGF=PGF,
               parm.names=parm.names,
               mon.names=mon.names,
#               pos.mu=pos.mu, pos.sigma=pos.sigma, pos.rho=pos.rho,
               N=N, yhat=ynew)

## prob2 <- function(parm,data){
##     ## mu
##     mu.prior <- sum(dnorm(parm[pmu],meanmu,stdmu,log=T))
##     ## sigma
##     parm[pls] <- interval(parm[pls],(lsigmaa),(lsigmab))
##     sigma.prior <- sum(dunif((parm[pls]),lsigmaa,lsigmab,log=T))
##     ## rho
##     U <- as(nearPD(as.symmetric.matrix(parm[pU],d))$mat,'matrix')
##     parm[pU] <-  lower.triangle(U,diag=T)
##     U.prior <- dinvwishart(U,nu,S,log=T)
##     ## likelihood
##     cho <- chol(U) %*% diag(exp(parm[pls])/sqrt(diag(U)))
##     LL <- sum(dmvnc(data$y, parm[pmu], cho, log=T))
##     ##
##     LP <- LL + mu.prior + sigma.prior + U.prior
##     return <- list(LP=LP, Dev=-2*LL, Monitor=1, yhat=1, parm=parm)
## }

prob <- function(parm,data){
    ## mu
    mu.prior <- sum(dnorm(parm[pmu],meanmu,stdmu,log=T))
    ## sigma
    parm[pls] <- interval(parm[pls],(lsigmaa),(lsigmab))
    sigma.prior <- sum(dunif((parm[pls]),lsigmaa,lsigmab,log=T))
    ## rho
    U <- as.symmetric.matrix(parm[pU],d)
    if(!is.positive.definite(U)) {
        U <- as.positive.definite(U)
        parm[pU] <-  lower.triangle(U,diag=T)
    }
    U.prior <- dinvwishart(U,nu,S,log=T)
    ## likelihood
    cho <- chol(U) %*% diag(exp(parm[pls])/sqrt(diag(U)))
    LL <- sum(dmvnc(data$y, parm[pmu], cho, log=T))
    ##
    LP <- LL + mu.prior + sigma.prior + U.prior
    return <- list(LP=LP, Dev=-2*LL, Monitor=1, yhat=1, parm=parm)
}

nchains <- 10
Initial.Values <- matrix(0,nchains,length(parm.names))
for(i in 1:nchains){
Initial.Values[i,] <- GIV(prob, mydata, n=1000, PGF=T)
}

## samplei <- LaplacesDemon(prob, mydata, GIV(prob, mydata, n=1000, PGF=T),
##                         Covar=NULL,
##                         Thinning=1,
##                         Iterations=10000, Status=1000,#LogFile=paste0(filename,'_LD_init_log'),
##                         Algorithm="AFSS", Specs=list(A=10000, B=NULL, m=100, n=0, w=1)
##                         )

sample2 <- LaplacesDemon.hpc(prob, mydata, Initial.Values,
                        Covar=NULL,
                        Thinning=10,
                        Iterations=1e6, Status=5e4,
                        Chains=nchains,CPUs=nchains,LogFile=paste0(filename,'_LDlog'), #Packages=c('Matrix'),#Type="MPI",
                        Algorithm="RDMH"#, Specs=list(B=list(1:d,d1:d2,d3:dnp))
                        ##Algorithm="Slice", Specs=list(B=list(1:d,d1:d2,d3:dnp), Bounds=list(c(-Inf,Inf), c(exp(lsigmaa),exp(lsigmab)), c(-500,500)), m=Inf, Type="Continuous", w=1)
                        ##Algorithm="AFSS", Specs=list(A=Inf, B=NULL, m=100, n=0, w=1)
)

sample0 <- Combine(sample2,mydata)

lml <- LML(Model=prob,Data=mydata,Modes=sample0$Summary1[1:dnp,'Median'],theta=sample0$Posterior1,LL=-sample0$Deviance/2,Covar=sample0$Covar,method='NSIS')[1]

nsamples <- dim(sample0$Posterior1)[1]
samples <- cbind(matrix(0,nsamples,d), sample0$Posterior1[,pmu], sample0$Posterior1[,pls],matrix(0,nsamples,nr))
for(i in 1:nsamples){
    U <- as.symmetric.matrix(sample0$Posterior1[i,pU],d)
    sdU <- sqrt(diag(U))
    cho <- chol(U)
    rhoc <- cho %*% diag(1/sdU)
    samples[i,(3*d+1):(3*d+nr)] <- upper.triangle(t(rhoc) %*% rhoc)
    ##
    samples[i,1:d] <- invlogit(rmvnc(1, sample0$Posterior1[i,pmu], cho %*% diag(exp(sample0$Posterior1[i,pls])/sdU)))
}

save.image(file=paste0(filename,'.RData'))


npar <- dim(samples)[2]
nbins <- 31
edg <- matrix(0,nbins+1,npar)
## xx <- matrix(0,npar,nbins)
## yy <- matrix(0,npar,nbins)
## yyc <- matrix(0,npar,nbins)

## functs <- list(
##     function(x){ifelse(x==0, 0, exp(-1/(2*x^2))/(x^3))},
##     function(x){ifelse(x==0, sqrt(2/pi),
##                        -1/(sqrt(2*pi)*(x^4)) +
##                        (exp(1/(2*x^2)) *
##                         (1+x^2)*
##                         (1- erf(1/(abs(x)*sqrt(2))))
##                        )/abs(2*x^5)
##                        )},
##     function(x){ifelse(x==0, 0, sqrt(2/pi)*exp(-1/(2*x^2))/(x^4))},
##     function(x){dunif(x,-1,1)}
## )

edg[,1:d] <- seq(0,1,length.out=nbins+1)
for(i in c(d+pmu, d+pls)){
edg[,i] <- seq(quantile(samples[,i],0.1), quantile(samples[,i],0.9) ,length.out=nbins+1)
}
##edg[,d+pls] <- seq((lsigmaa),(lsigmab),length.out=nbins+1)##seq(min(rsamples[,i]), max(rsamples[,i]), length.out=nbins+1)
edg[,(3*d)+(1:nr)] <- seq(-1,1,length.out=nbins+1)##seq(min(rsamples[,i]), max(rsamples[,i]), length.out=nbins+1)

## for(i in 1:npar){
## xx[i,] <- (edg[-1,i] + edg[-length(edg[i,,i])])/2
## yy[i,] <- functs[[i]](xx[i,])
## yyc[i,] <- yy[i,]*nsamples*(edg[2,i]-edg[1,i])
## }


pdf(paste0(filename,'_dens.pdf'))
plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
text(1, 9, ptitle, pos=4,cex=0.5)
text(1, 8, paste0(d,' graph prop, ',N,' indiv'), pos=4,cex=0.5)
text(1, 7, paste0('LML = ',lml), pos=4,cex=0.5)
#
for(i in 1:d){
    posterior <- samples[,i]
    dat <- data.frame(x=posterior)
   # data.f <- data.frame(x=xx[i,], y=yy[i,])
    print(
        ggplot() + geom_histogram(data=dat,aes(x=x,y=..density..),
                                  breaks=edg[,i]) #+ scale_y_log10()
        +labs(x=tparm.names[i],title=paste0('mean =',signif(mean(posterior),3),', sd = ',signif(sd(posterior),3)))
        + geom_vline(xintercept=mean(posterior),colour=mygreen)
        + geom_vline(xintercept=median(posterior),colour=mygreen,linetype='dashed')
        #+ geom_vline(xintercept=datam[,i],colour=myred)
        + geom_point(aes(x=datam[,i], y=0),size=2,colour=myred,shape=16)
        #+xlim(-3*sd(posterior),3*sd(posterior))
    )
}
for(i in d1:npar){
    posterior <- samples[,i]
    dat <- data.frame(x=posterior)
   # data.f <- data.frame(x=xx[i,], y=yy[i,])
    print(
        ggplot() + geom_histogram(data=dat,aes(x=x,y=..density..),
                                  breaks=edg[,i]) #+ scale_y_log10()
        +labs(x=tparm.names[i],title=paste0('mean =',signif(mean(posterior),3),', sd = ',signif(sd(posterior),3)))
        + geom_vline(xintercept=mean(posterior),colour=mygreen)
        + geom_vline(xintercept=median(posterior),colour=mygreen,linetype='dashed')
        #+xlim(-3*sd(posterior),3*sd(posterior))
    )
}
for(j in 1:(d-1)){
    for(i in (j+1):d){
        posterior <- samples[,j]
        posterior2 <- samples[,i]
        dat <- data.frame(x=posterior,y=posterior2)
        gmedian <- Gmedian(cbind(posterior,posterior2))
        datapoints <- data.frame(x=datam[,j], y=datam[,i])
        print(
            ggplot() + stat_bin2d(data=dat,aes(x=x,y=y,fill=..density..),
                                  breaks=list(x=edg[,j], y=edg[,i])
                                  ) +
            scale_fill_gradientn(colours=brewer.pal(n=9,name='Blues')) +
            labs(x=tparm.names[j],y=tparm.names[i],title=ptitle) 
            +geom_point(data=datapoints ,aes(x,y),colour=myred,size=4,shape=16)
            + geom_point(aes(x=mean(posterior), y=mean(posterior2)),size=3,colour=mygreen,shape=17)
            + geom_point(aes(x=gmedian[1], y=gmedian[2]),size=3,colour=mygreen,shape=15)
            + coord_fixed()
            #+xlim(-3*sd(posterior),3*sd(posterior))
            #+ylim(-3*sd(posterior2),3*sd(posterior2))
        )
        }
}
dev.off()    
##logs
pdf(paste0(filename,'_logcounts.pdf'))
plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
text(1, 9, ptitle, pos=4,cex=0.5)
text(1, 8, paste0(d,' graph prop, ',N,' indiv'), pos=4,cex=0.5)
text(1, 7, paste0('LML = ',lml), pos=4,cex=0.5)
#
for(i in 1:d){
    posterior <- samples[,i]
    dat <- data.frame(x=posterior)
    #data.f <- data.frame(x=xx[i,], y=yyc[i,])
    print(
        ggplot() + geom_histogram(data=dat,aes(x=x,y=..count..),
                                  breaks=edg[,i])+ scale_y_log10()
        +labs(x=tparm.names[i],title=paste0('mean =',signif(mean(posterior),3),', sd = ',signif(sd(posterior),3)))
        + geom_vline(xintercept=mean(posterior),colour=myyellow)
        + geom_vline(xintercept=median(posterior),colour=myyellow,linetype='dashed')
        #+ geom_vline(xintercept=datam[,i],colour=myred)
        + geom_point(aes(x=datam[,i], y=1),size=2,colour=myred,shape=16)
        #+xlim(-3*sd(posterior),3*sd(posterior))
    )
}
for(i in d1:npar){
    posterior <- samples[,i]
    dat <- data.frame(x=posterior)
    #data.f <- data.frame(x=xx[i,], y=yyc[i,])
    print(
        ggplot() + geom_histogram(data=dat,aes(x=x,y=..count..),
                                  breaks=edg[,i])+ scale_y_log10()
        +labs(x=tparm.names[i],title=paste0('mean =',signif(mean(posterior),3),', sd = ',signif(sd(posterior),3)))
        + geom_vline(xintercept=mean(posterior),colour=myyellow)
        + geom_vline(xintercept=median(posterior),colour=myyellow,linetype='dashed')
        #+xlim(-3*sd(posterior),3*sd(posterior))
    )
}
for(j in 1:(d-1)){
    for(i in (j+1):d){
        posterior <- samples[,j]
        posterior2 <- samples[,i]
        dat <- data.frame(x=posterior,y=posterior2)
        datapoints <- data.frame(x=datam[,j], y=datam[,i])
        gmedian <- Gmedian(cbind(posterior,posterior2))
        print(
            ggplot() + stat_bin2d(data=dat,aes(x=x,y=y,fill=..count..),
                                  breaks=list(x=edg[,j], y=edg[,i])
                                  ) +
            scale_fill_gradientn(colours=brewer.pal(n=9,name='Blues'),trans='log10') +
            labs(x=tparm.names[j],y=tparm.names[i],title=ptitle) 
            +geom_point(data=datapoints ,aes(x,y),colour=myred,size=4,shape=16)
            + geom_point(aes(x=mean(posterior), y=mean(posterior2)),size=3,colour=mygreen,shape=17)
            + geom_point(aes(x=gmedian[1], y=gmedian[2]),size=3,colour=mygreen,shape=15)
            + coord_fixed()
            #+xlim(-3*sd(posterior),3*sd(posterior))
            #+ylim(-3*sd(posterior2),3*sd(posterior2))

        )
        }
    }
dev.off()

plot(sample0, BurnIn=0, mydata, PDF=T, Parms=NULL)
file.rename('LaplacesDemon.Plots.pdf',paste0('LDPlots_',filename,'.pdf'))
## png(paste0('LDPlots_',filename,'.png'))
## plot(sample0, BurnIn=0, mydata, PDF=F, Parms=NULL)
## dev.off()

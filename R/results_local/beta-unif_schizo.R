## test script for normal model 
library('MCMCpack')
library('LaplacesDemon')
library('RColorBrewer')
library('bayesplot')
library('ggplot2')
library('Matrix')
library('Gmedian')
library('rootSolve')
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
filename <- 'beta-unif_schizo_all_rdmh'
ptitle <- 'beta-uniform model, Schizo all (RDMH)'
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
dnp <- 2*d # total num. parameters
d1 <- d+1
d2 <- 2*d
## suff statistics: log-means and log-1-means
lmean <- colMeans(log(datam))
l1mean <- colMeans(log(1-datam))
sdlm <- sd(lmean)
sdl1m <- sd(l1mean)
fip <- function(x,y){return(digamma(x)-y-digamma(sum(x)))}
initparms <- matrix(0,2,d)
for(i in 1:d){
    initparms[,i] <- multiroot(f=fip,start=rep(1,2),positive=T,parms=c(lmean[i],l1mean[i]))$root
}

## new datum - we take the mean of the data we have just for check
ynew <- dmean

## hyperparameters for hyperprior
meana <- 0 # mean for log-a
loa <- -5
upa <- 5
meanb <- 0 # mean for log-b
lob <- -5
upb <- 5

## parameters:  means + elements of Cholesky decomposition
pa <- 1:d
pb <- d1:d2
parm.names <- as.parm.names(list(a=rep(0,d), b=rep(0,d)))

# parameters we want to see
tparm.names <- c(y.names, parm.names)

mon.names <- c('')

## function to generate initial values
## PGF2 <- function(data){
##     mu <- rnorm(d,meanmu,stdmu)
##     lsigma <- runif(d,lsigmaa,lsigmab)
##     U <- rinvwishart(d+1,diag(d))
##     return(c(mu, lsigma, lower.triangle(U,diag=T)))
## }
PGF <- function(data){
    return(c(
    interval(rnorm(d, log(initparms[1,]),3*sd(log(initparms[1,]))),loa,upa),
    interval(rnorm(d, log(initparms[2,]),3*sd(log(initparms[2,]))),loa,upa)
    ))
}

mydata <- list(y=datam, PGF=PGF,
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
    parm <- interval(parm,loa,upa)
    ab.prior <- sum(dunif(parm,loa,upa,log=T))
    LL <- sum(dbeta(t(data$y), exp(parm[pa]), exp(parm[pb]), log=T))
    ##
    LP <- LL + ab.prior
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
                        Thinning=1,
                        Iterations=50000, Status=1000,
                        Chains=nchains,CPUs=nchains,LogFile=paste0(filename,'_LDlog'), #Packages=c('Matrix'),#Type="MPI",
                        Algorithm="RDMH"#, Specs=list(B=list(1:d,d1:d2,d3:dnp))
                        ##Algorithm="Slice", Specs=list(B=list(1:d,d1:d2,d3:dnp), Bounds=list(c(-Inf,Inf), c(exp(lsigmaa),exp(lsigmab)), c(-500,500)), m=Inf, Type="Continuous", w=1)
                        ##Algorithm="AFSS", Specs=list(A=Inf, B=NULL, m=100, n=0, w=1)
)

sample0 <- Combine(sample2,mydata)

lml <- LML(Model=prob,Data=mydata,Modes=sample0$Summary1[1:dnp,'Median'],theta=sample0$Posterior1,LL=-sample0$Deviance/2,Covar=sample0$Covar,method='NSIS')[1]

nsamples <- dim(sample0$Posterior1)[1]
samples <- cbind(matrix(0,nsamples,d), sample0$Posterior1)
for(i in 1:nsamples){
    samples[i,1:d] <- rbeta(d,exp(sample0$Posterior1[i,pa]),exp(sample0$Posterior1[i,pb]))
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
for(i in (d+1):(d+dnp)){
edg[,i] <- seq(quantile(samples[,i],0.1), quantile(samples[,i],0.9) ,length.out=nbins+1)
}

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

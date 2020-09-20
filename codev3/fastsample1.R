## Monte Carlo calculation of posterior probability ##

ntasks <- as.numeric(commandArgs(trailingOnly=TRUE))

pdff <- function(filename){pdf(file=paste0(filename,'.pdf'),paper='a4r',height=11.7,width=16.5)}
#library('ggplot2')
library('RColorBrewer')
#library('cowplot')
#library('png')
#library('plot3D')
library('foreach')
library('doFuture')

registerDoFuture()
print(paste0('available workers: ', availableCores()))
if(length(ntasks)==0){ntasks <- availableCores()}
plan(multiprocess, workers=ntasks-1)
print(paste0('number of workers: ', nbrOfWorkers()))

library('doRNG')
#library('LaplacesDemon')
library('dirichletprocess')
library('mvtnorm')
library('RNetCDF')
#library('Rmpfr')
options(bitmapType='cairo')
mypurpleblue <- '#4477AA'
myblue <- '#66CCEE'
mygreen <- '#228833'
myyellow <- '#CCBB44'
myred <- '#EE6677'
myredpurple <- '#AA3377'
mygrey <- '#BBBBBB'
mycolours <- c(myblue, myred, mygreen, myyellow, myredpurple, mypurpleblue, mygrey, 'black')
palette(mycolours)
barpalette <- colorRampPalette(c(mypurpleblue,'white',myredpurple),space='Lab')
barpalettepos <- colorRampPalette(c('white','black'),space='Lab')
dev.off()

set.seed(181225+0)

#################################################################
######################## EQUATIONS SETUP ########################
#################################################################

generalname <- paste0('testDP')
datafilename <- 'data1'

nData <- 10;
data <- matrix(rnorm(2*nData,mean=c(-1,-1,1,1),sd=0.1), nrow=nData, ncol=2, byrow=TRUE)

dp <- DirichletProcessMvnormal(data)

nSamples <- 1000
fitDp <- Fit(dp, its=nSamples)



## p(ynew[m] | data, hyperp)
postmargx <- function(x,dpobj,m){
    if(is.null(dim(x))){dim(x) <- c(length(x),1)}
    n <- dpobj$n
    
    G <- dpobj$mixingDistribution$priorParameters
    df <- G$nu - length(G$mu0) + 1
    sigma <- G$Lambda * (G$kappa0 + 1)/(G$kappa0 * df)


    ac <- dpobj$alphaChain
    wc <- dpobj$weightsChain
    parsc <- dpobj$clusterParametersChain
    nsamples <- length(wc)
    
    foreach(i=seq_len(nsamples), .combine='+')%dopar%{
        w <- wc[[i]]
        pars <- parsc[[i]]
        foreach(k=seq_len(length(w)), .combine='+')%do%{
            dnorm(x=x, mean=pars$mu[1,m,k], sd=sqrt(pars$sig[m,m,k]), log=FALSE) * w[k]
        } /(ac[i] + n)
    } * n/nsamples + mean(ac/(ac+n)) * dmvt(x=x, delta=c(G$mu0[m]), sigma=matrix(sigma[m,m],1,1), df=df, type='shifted', log=FALSE)
}

## p(ynew[m] | data, hyperp) WITHOUT t-term (for checks)
postdmargx <- function(x,dpobj,m){
    if(is.null(dim(x))){dim(x) <- c(length(x),1)}
    n <- dpobj$n
    
    G <- dpobj$mixingDistribution$priorParameters
    df <- G$nu - length(G$mu0) + 1
    sigma <- G$Lambda * (G$kappa0 + 1)/(G$kappa0 * df)


    ac <- dpobj$alphaChain
    wc <- dpobj$weightsChain
    parsc <- dpobj$clusterParametersChain
    nsamples <- length(wc)
    
    foreach(i=seq_len(nsamples), .combine='+')%dopar%{
        w <- wc[[i]]
        pars <- parsc[[i]]
        foreach(k=seq_len(length(w)), .combine='+')%do%{
            w[k] * dnorm(x=x, mean=pars$mu[1,m,k], sd=sqrt(pars$sig[m,m,k]), log=FALSE)
        } 
    } /nsamples
}

## p(ynew | data, hyperp) - fastest version
posteriorx <- function(x,dpobj){
    if(is.null(dim(x))){dim(x) <- c(1,length(x))}
    n <- dpobj$n
    
    G <- dpobj$mixingDistribution$priorParameters
    df <- G$nu - length(G$mu0) + 1
    sigma <- G$Lambda * (G$kappa0 + 1)/(G$kappa0 * df)

    ac <- dpobj$alphaChain
    nsamples <- length(ac)
    
   c(foreach(i=seq_len(nsamples), .combine='+')%dopar%{
        LikelihoodFunction(dpobj=dpobj, ind=i)(x)/(ac[i] + n)
        }) * n/nsamples + mean(ac/(ac+n)) * dmvt(x=x, delta=G$mu0, sigma=sigma, df=df, type='shifted', log=FALSE)
}

## p(ynew | data, hyperp) - built-in version
rpostx <- function(x,dpobj){
    nsamples <- length(dpobj$alphaChain)
    if(is.null(dim(x))){dim(x) <- c(1,length(x))}
    
    c(foreach(i=seq_len(nsamples), .combine='+')%dopar%{
        PosteriorFunction(dpobj, ind=i)(x)
    })/nsamples
}

llf <- function(x,dpobj,ind){
    n <- dpobj$n
    
    G <- dpobj$mixingDistribution$priorParameters
    df <- G$nu - length(G$mu0) + 1
    sigma <- G$Lambda * (G$kappa0 + 1)/(G$kappa0 * df)


    ac <- dpobj$alphaChain
    wc <- dpobj$weightsChain
    parsc <- dpobj$clusterParametersChain
    nsamples <- length(wc)
    
    w <- wc[[ind]]
    pars <- parsc[[ind]]
    foreach(j=seq_len(length(w)), .combine='+')%do%{
        dmvnorm(x=x, mean=pars$mu[,,j], sigma=pars$sig[,,j], checkSymmetry=FALSE, log=FALSE) * w[j]
    }
}


## test MC chain
testdp <- dp
ndata <- testdp$n
nsamples <- 5
set.seed(111)
tests <- foreach(i=1:nsamples)%do%{
    testdp <- UpdateAlpha(ClusterParameterUpdate(ClusterComponentUpdate(testdp)))
    list(a=testdp$alpha, w=tabulate(testdp$clusterLabels)/ndata, t=testdp$clusterParameters)
}

## p(ynew[m] | data, hyperp) - via sampling from base measure
pmarg <- function(x,dpobj,m,ind=length(dpobj$alphaChain)){clu <- PosteriorClusters(dpobj,ind)
    foreach(i=1:length(clu$weights), .combine='+')%dopar%{
    clu$weights[i] * dnorm(x=x, mean=clu$params[[1]][1,m,i], sd=sqrt(clu$params[[2]][m,m,i]), log=FALSE)
}}


postmargxi <- function(x,dpobj,m,i=length(dpobj$alphaChain)){
    if(is.null(dim(x))){dim(x) <- c(length(x),1)}
    n <- dpobj$n
    
    G <- dpobj$mixingDistribution$priorParameters
    df <- G$nu - length(G$mu0) + 1
    sigma <- G$Lambda * (G$kappa0 + 1)/(G$kappa0 * df)


    ac <- dpobj$alphaChain
    wc <- dpobj$weightsChain
    parsc <- dpobj$clusterParametersChain
    nsamples <- length(wc)
    
        w <- wc[[i]]
        pars <- parsc[[i]]
        foreach(k=seq_len(length(w)), .combine='+')%do%{
        w[k] * dnorm(x=x, mean=pars$mu[1,m,k], sd=sqrt(pars$sig[m,m,k]), log=FALSE)
        } *n/(ac[i] + n) +
            ac[i]/(ac[i] + n) * dmvt(x=x, delta=c(G$mu0[m]), sigma=sigma[m,m]*diag(1), df=df, type='shifted', log=FALSE)
}



xfitdp <- fitdp2
nsamples <- length(xfitdp$alphaChain)
sseq <- round(seq(1,nsamples,length.out=100))
xg <- seq(-5,5,by=0.01)
matplot(xg,postmargxi(xg,xfitdp,m=1,sseq[1]),type='l',ylim=c(0,0.5))
for(i in sseq[-1]){
matplot(xg,postmargxi(xg,xfitdp,m=1,i),type='l', add=TRUE)
}

xfitdp <- fitdp2
xg <- seq(-5,5,by=0.01)
matplot(xg,postmargx(xg,xfitdp,1),type='l',ylim=c(0,0.5))



logit2 <- function(x){log(1+x)-log(1-x)}
jlogit <- function(x){2/((1+x)*(1-x))}

datasize <- 10
y <- rbind(rmvnorm(n=datasize, mean=c(-2,-2), sigma=0.5^2 * diag(2)),
           rmvnorm(n=datasize, mean=c(2,2), sigma=0.5^2 * diag(2)))


## healthy
datafilename <- 'weights_con_40cons'
y <- t(as.matrix(read.csv(paste0(datafilename,'.csv'),header=FALSE,sep=',')))
dims <- ncol(y)
y <- logit2(y)

dp <- DirichletProcessMvnormal(y,
                                #alphaPriors=c(0.1, 0.1),
                                g0Priors=list(mu0=rep(0,dims), kappa0=0.01, Lambda=diag(dims)/2, nu=dims+1), numInitialClusters = 5)
fitdp <- dp

for(i in 1:2){fitdp <- Fit(fitdp,5000)}


## Schizo
datafilename <- 'weights_Schizo_40cons'
y2 <- t(as.matrix(read.csv(paste0(datafilename,'.csv'),header=FALSE,sep=',')))
dims <- ncol(y2)
y2 <- logit2(y2)


dp2 <- DirichletProcessMvnormal(y2,
                                #alphaPriors=c(0.1, 0.1),
                                g0Priors=list(mu0=rep(0,dims), kappa0=0.01, Lambda=diag(dims)/2, nu=dims+1), numInitialClusters = 5)
fitdp2 <- dp2

for(i in 1:2){fitdp2 <- Fit(fitdp2,5000)}


## plots
for(marg in 1:dims){
    pdff(paste0('H_vs_S_marg',marg))
xfitdp <- fitdp
nsamples <- length(xfitdp$alphaChain)
sseq <- round(seq(1,nsamples,length.out=50))
xg <- seq(-1,1,by=0.01)
matplot(xg,jlogit(xg)*pmarg(logit2(xg),xfitdp,m=marg,sseq[1]),type='l',col=myblue,ylim=c(0,1.5), xlab=paste0('correl. # ',marg), ylab='prob density')
for(i in sseq[-1]){
matplot(xg,jlogit(xg)*pmarg(logit2(xg),xfitdp,m=marg,i),type='l',col=myblue, add=TRUE)
}
xfitdp <- fitdp2
nsamples <- length(xfitdp$alphaChain)
sseq <- round(seq(1,nsamples,length.out=50))
matplot(xg,jlogit(xg)*pmarg(logit2(xg),xfitdp,m=marg,sseq[1]),type='l',col=myred,ylim=c(0,1.5), add=TRUE)
for(i in sseq[-1]){
matplot(xg,jlogit(xg)*pmarg(logit2(xg),xfitdp,m=marg,i),type='l',col=myred, add=TRUE)
}
    dev.off()
}




##############################################################################
##############################################################################
##############################################################################


## read example sample
dataabsfreqs <- as.matrix(read.csv(paste0(datafilename,'.csv'),header=FALSE,sep=','))
## format: 1 row = 1 stimulus, 1 column = 1 response

nmcsamples <- 1e5

nstimuli <- nrow(dataabsfreqs)
nresponses <- ncol(dataabsfreqs)
datasize <- rowSums(dataabsfreqs)
datarelfreqs <- dataabsfreqs/datasize

## Mutual info: argument is matrix of *relative* conditional frequencies:
## 1 row = 1 stimulus, 1 column = 1 response
## output is in base=number of stimuli (2 stimuli -> bits)
MI <- function(condfreqs){
    nstim <- nrow(condfreqs)
    sum(condfreqs*log(t(t(condfreqs)/colMeans(condfreqs)), base=nstim), na.rm=TRUE)/nstim
}

## MI of the example sample
dataMI <- MI(datarelfreqs)

## Dirichlet parameters for F-uniform prior
priorsize <- nresponses * 10

dirichpriorfreqs <- foreach(i=1:nstimuli, .combine=rbind)%do%{rep(1/nresponses, nresponses)}
dirichpriorsize <- foreach(i=1:nstimuli, .combine=c)%do%{priorsize}

dirichpostparams <- dirichpriorsize * dirichpriorfreqs + dataabsfreqs

## generate posterior samples of long-run frequency distrs
Fsamples <- aperm(simplify2array(foreach(i=1:nstimuli)%do%{
    rdirichlet(n=nmcsamples, alpha=dirichpostparams[i,]) }) , perm=c(1,3,2))
## format: 1 row = 1 MCsample, 1 column = 1 stimulus, 1 3rdDim = 1 response

## generate posterior samples of corresponding long-run MI
MIsamples <- foreach(i=1:nmcsamples, .combine=c)%dopar%{MI(Fsamples[i,,])}

## generate posterior samples of corresponding log-likelihoods
#LLLsamples <- foreach(i=1:nmcsamples, .combine=c)%:%foreach(s=1:nstimuli, .combine='+')%dopar%{dmultinom(dataabsfreqs[s,], prob=Fsamples[i,s,], log=TRUE)}
LLLsamples <- foreach(i=1:nmcsamples, .combine=c)%:%foreach(s=1:nstimuli, .combine='+')%dopar%{lfactorial(datasize[s]) + sum(dataabsfreqs[s,]*log(Fsamples[i,s,])-lfactorial(dataabsfreqs[s,]))}
    
## save posterior samples
#saveRDS(Fsamples, paste0('Fsamples_', generalname, priorsize,'.rds'))
saveRDS(MIsamples, paste0('MIsamples_', generalname, priorsize,'.rds'))
saveRDS(LLLsamples, paste0('LLLsamples_', generalname, priorsize,'.rds'))

write.table(cbind(LLLsamples, MIsamples),file= paste0('LLL-MI_samples_', generalname, priorsize,'.csv'),row.names=FALSE,col.names=FALSE, sep=',')

## find and same long-run frequencies with min & max MI
#minmi <- which(MIsamples==min(MIsamples))
#maxmi <- which(MIsamples==max(MIsamples))

## write.table(Fsamples[minmi,,],file= paste0('minMI_Fsamples_', generalname, priorsize,'.csv'),row.names=FALSE,col.names=FALSE, sep=',')
## write.table(Fsamples[maxmi,,],file= paste0('maxMI_Fsamples_', generalname, priorsize,'.csv'),row.names=FALSE,col.names=FALSE, sep=',')


exit()

## credibility interval & median
quants <- unname(quantile(MIsamples, probs=c(0.025, 0.5, 0.975)))

## histogram of posterior for long-run MI
histcells <- seq(from=0, to=1, length.out=round(75/diff(range(MIsamples))))
pdff(paste0(generalname,priorsize))
hist(MIsamples, breaks=histcells, freq=FALSE, xlab='long-run MI/bit',
     main=paste0('Dirichlet prior with uniform reference distribution and prior size = ', priorsize, '\n95% probability range: (', signif(quants[1],3), ', ', signif(quants[3],3), ') bit\nmedian: ', signif(quants[2],3),
                 ' bit\nsample MI = ', signif(dataMI,3), ' bit'))
matpoints(x=dataMI, y=0, type='p', pch=17, cex=2, col=myred)
dev.off()

#####
## Forcasts using only the mutual info of the sample, rather than the sample frequencies

nmcsamplesbis <- 40000

## Dirichlet parameters for F-uniform prior
dirichpriorfreqs <- foreach(i=1:nstimuli, .combine=rbind)%do%{rep(1/nresponses, nresponses)}
dirichpriorsize <- foreach(i=1:nstimuli, .combine=c)%do%{nresponses}*10

dirichpriorparams <- dirichpriorsize * dirichpriorfreqs

Fpriorsamples <- aperm(simplify2array(foreach(i=1:nstimuli)%do%{
    rdirichlet(n=nmcsamplesbis, alpha=dirichpriorparams[i,]) }) , perm=c(1,3,2))
## format: 1 row = 1 MCsample, 1 column = 1 stimulus, 1 3rdDim = 1 response

doubleMIsamples <- foreach(i=1:nmcsamplesbis, .combine=rbind)%dorng%{
    c(MI(foreach(s=1:nstimuli, .combine=rbind)%do%{
        tabulate(sample(x=1:nresponses, size=datasize[s], replace=TRUE,
                        Fpriorsamples[i,s,]), nbins=nresponses)/datasize[s]
    })
    , MI(Fpriorsamples[i,,]))}


miwidth <- 0.01 # width around the mutual info
MIsamplesbis <- doubleMIsamples[(doubleMIsamples[,1] > dataMI - miwidth) & (doubleMIsamples[,1] < dataMI + miwidth),2]

histcells <- seq(from=0, to=1, length.out=round(20/diff(range(MIsamplesbis))))
pdff(paste0('onlyMI_',generalname,priorsize))
hist(MIsamplesbis, breaks=histcells, freq=FALSE, xlab='long-run MI/bit',
     main=paste0('Dirichlet prior with uniform reference distribution and prior size = ', priorsize, '\nsample MI = ', signif(dataMI,3), ' bit'))
matpoints(x=dataMI, y=0, type='p', pch=17, cex=2, col=myred)
dev.off()


mcdataMI <- doubleMIsamples[,1]
histcells <- seq(from=0, to=1, length.out=round(25/diff(range(mcdataMI))))
pdff(paste0('test',priorsize))
hist(mcdataMI, breaks=histcells, freq=FALSE, xlab='long-run MI/bit',
     main=paste0('Dirichlet prior with uniform reference distribution and prior size = ', priorsize, '\nsample MI = ', signif(dataMI,3), ' bit'))
matpoints(x=dataMI, y=0, type='p', pch=17, cex=2, col=myred)
dev.off()


plan(sequential)

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
set.seed(181225+111)

#################################################################
######################## EQUATIONS SETUP ########################
#################################################################
session <- '001001'

groups <- list(H='con', S='sch')
numH <- 55
numS <- 49

HfromH <- foreach(s = 1:numH, .combine=c)%do%{
    unname(unlist(read.csv(paste0('logscore',session,'_', groups$H, '_marg_', s, '.csv'),header=FALSE,sep=',')))
}

HfromS <- unname(unlist(read.csv(paste0('logscore',session,'_', groups$S, '_marg_', numS+1, '.csv'),header=FALSE,sep=',')))

looH <- cbind(100*exp(HfromH)/(exp(HfromH)+exp(HfromS)), HfromH - HfromS, HfromH, HfromS)

write.table(looH, file= paste0('cvls_', groups$H,'.csv'), row.names=FALSE, col.names=FALSE, sep=',')


SfromS <- foreach(s = 1:numS, .combine=c)%do%{
    unname(unlist(read.csv(paste0('logscore',session,'_', groups$S, '_marg_', s, '.csv'),header=FALSE,sep=',')))
}

SfromH <- unname(unlist(read.csv(paste0('logscore',session,'_', groups$H, '_marg_', numH+1, '.csv'),header=FALSE,sep=',')))

looS <- cbind(100*exp(SfromS)/(exp(SfromS)+exp(SfromH)), SfromS - SfromH, SfromS, SfromH)

write.table(looS, file= paste0('cvls_', groups$S,'.csv'), row.names=FALSE, col.names=FALSE, sep=',')


filesink <- file('cvls_means.txt')
sink(filesink)
print(paste0('mean cvls H: ', colMeans(looH)))
print(paste0('mean cvls S: ', colMeans(looS)))

print(paste0('mean cvls all: ', colMeans(rbind(looH, looS))))

print(paste0('hits cvls H: ', colMeans(looH > 0)))
print(paste0('hits cvls S: ', colMeans(looS > 0)))

print(paste0('hits cvls all: ', colMeans(rbind(looH, looS) > 0)))
sink()




hist(looS[,1])







 <- paste0('fits_', groups, '_marg_')
datafilename <- paste0('weights_', group, '_40cons')

yH <- t(as.matrix(read.csv(paste0('weights_', groups$H, '_40cons.csv'),header=FALSE,sep=',')))
datasizeH <- nrow(yH)
yH <- logit2(yH)

yS <- t(as.matrix(read.csv(paste0('weights_', groups$S, '_40cons.csv'),header=FALSE,sep=',')))
datasizeS <- nrow(yS)
yS <- logit2(yS)

dims <- ncol(yH)

## p(ynew | data, hyperp) - fastest version
posteriorx <- function(x,dpobj){
    if(is.null(dim(x))){dim(x) <- c(1,length(x))}
    n <- dpobj$n
    
    G <- dpobj$mixingDistribution$priorParameters
    df <- G$nu - length(G$mu0) + 1
    sigma <- G$Lambda * (G$kappa0 + 1)/(G$kappa0 * df)

    ac <- dpobj$alphaChain
    nsamples <- length(ac)
    
   c(foreach(i=seq_len(nsamples), .combine='+')%do%{
        LikelihoodFunction(dpobj=dpobj, ind=i)(x)/(ac[i] + n)
        }) * n/nsamples + mean(ac/(ac+n)) * dmvt(x=x, delta=G$mu0, sigma=sigma, df=df, type='shifted', log=FALSE)
}


system.time(for(i in 1:1){fitdpS <- readRDS(paste0('fits_', groups$S, '_marg_', datasizeS+1,'.rds')) ;pxH <- posteriorx(yH,fitdpS)})

system.time(for(i in 1:1){
looH <- foreach(s = 1:1, .combine=rbind)%do%{
    x <- yH[s,]
    dim(x) <- c(1,dims)

    fitdpH <- readRDS(paste0('fits_', groups$H, '_marg_', s,'.rds'))

    log(c(posteriorx(x,fitdpH), pxH[s]))
}
})

write.table(cbind(looH[,1] - looH[,2], looH), file= paste0('cvls_', group[1],'.csv'), row.names=FALSE, col.names=FALSE, sep=',')


fitdpH <- readRDS(paste0('fits_', groups$H, '_marg_', datasizeH+1,'.rds'))
pxS <- posteriorx(yS,fitdpH)
looS <- foreach(s = 1:datasizeS, .combine=rbind)%do%{
    x <- yS[s,]
    dim(x) <- c(1,dims)

    fitdpS <- readRDS(paste0('fits_', groups$S, '_marg_', s,'.rds'))

    log(c(posteriorx(x,fitdpS), pxS[s]))
}

write.table(cbind(looS[,1] - looS[,2], looS), file= paste0('cvls_', group[2],'.csv'), row.names=FALSE, col.names=FALSE, sep=',')

pdff('cvls_means')
print(paste0('mean cvls H: ' mean(looH)))
print(paste0('mean cvls S: ' mean(looS)))

print(paste0('mean cvls all: ' mean(c(looH, looS)))

print(paste0('hits cvls H: ' mean(looH > 0)))
print(paste0('hits cvls S: ' mean(looS > 0)))

print(paste0('hits cvls all: ' mean(c(looH, looS) > 0))
dev.off()



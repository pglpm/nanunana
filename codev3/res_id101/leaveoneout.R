## Monte Carlo calculation of posterior probability ##

#ntasks <- as.numeric(commandArgs(trailingOnly=TRUE))

pdff <- function(filename){pdf(file=paste0(filename,'.pdf'),paper='a4r',height=11.7,width=16.5)}
#library('ggplot2')
library('RColorBrewer')
#library('cowplot')
#library('png')
#library('plot3D')
library('foreach')
library('doFuture')

## registerDoFuture()
## print(paste0('available workers: ', availableCores()))
## if(length(ntasks)==0){ntasks <- availableCores()}
## plan(multiprocess, workers=ntasks-1)
## print(paste0('number of workers: ', nbrOfWorkers()))

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
#dev.off()
set.seed(181225+111)

#################################################################
######################## EQUATIONS SETUP ########################
#################################################################
basename <- 'id101'

groups <- list(H='con', S='sch')
numH <- 55
numS <- 49

HfromH <- foreach(s = 1:numH, .combine=c)%do%{
    unname(unlist(read.csv(paste0(basename,'_logscore_', groups$H, '_', s, '.csv'),header=FALSE,sep=',')))
}

HfromS <- unname(unlist(read.csv(paste0(basename,'_logscore_', groups$S, '_', numS+1, '.csv'),header=FALSE,sep=',')))

looH <- cbind(100*exp(HfromH)/(exp(HfromH)+exp(HfromS)), HfromH - HfromS, HfromH, HfromS)

write.table(looH, file= paste0('cvls_', groups$H,'.csv'), row.names=FALSE, col.names=FALSE, sep=',')


SfromS <- foreach(s = 1:numS, .combine=c)%do%{
    unname(unlist(read.csv(paste0(basename,'_logscore_', groups$S, '_', s, '.csv'),header=FALSE,sep=',')))
}

SfromH <- unname(unlist(read.csv(paste0(basename,'_logscore_', groups$H, '_', numH+1, '.csv'),header=FALSE,sep=',')))

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

pdff(paste0('p_histogs'))
hist(looH[,1],breaks=seq(0,110,by=10)-5)
hist(looS[,1],breaks=seq(0,110,by=10)-5)
dev.off()

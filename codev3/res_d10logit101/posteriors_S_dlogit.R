## Monte Carlo calculation of posterior probability ##

s <- as.numeric(commandArgs(trailingOnly=TRUE))
ntasks <- 21

pdff <- function(filename){pdf(file=paste0(filename,'.pdf'),paper='a4r',height=11.7,width=16.5)}
#library('ggplot2')
library('RColorBrewer')
#library('cowplot')
#library('png')
#library('plot3D')
library('foreach')
library('doFuture')
## registerDoFuture()
## #library("future.batchtools")

## print(paste0('available workers: ', availableCores()))
## if(length(ntasks)==0){ntasks <- availableCores()}
## #plan(batchtools_slurm, workers=ntasks-1)
## #plan(multiprocess, workers=ntasks-1)
## cl <- makeCluster(ntasks-1)
## plan(cluster, workers = cl)
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
dev.off()

set.seed(110+s)

#################################################################
######################## EQUATIONS SETUP ########################
#################################################################
source('fdefslogit.R')
basename <- 'd10logit101'

rdims <- 1:10

groups <- c('con', 'sch')
group <- groups[2]

datafilename <- paste0('weights_', group, '_40cons')

y <- t(as.matrix(read.csv(paste0(datafilename,'.csv'),header=FALSE,sep=',')))
y <- y[,rdims]
dims <- ncol(y)
y <- logit2(y)
datasize <- nrow(y)

#foreach(s = 1:(datasize+1))%dopar%{
    print(paste0('marg: ',s))
    x <- y[-s,]
    
    dp <- DirichletProcessMvnormal(x,
                                        alphaPriors=c(1, 0.1),
                               g0Priors=list(mu0=rep(0,dims), kappa0=1, Lambda=diag(dims)/(dims+2), nu=dims+1), numInitialClusters = 3)

savefile <- paste0(basename,'_fits_',group,'_', s, '.rds')
if(file.exists(savefile)){
    fitdp <- readRDS(savefile)
} else {fitdp <- dp}
#fitdp <- myburn(fitdp,its=10000,progressBar=FALSE)

fitdp <- myfitfast(fitdp,its=1000,progressBar=FALSE)
saveRDS(fitdp,savefile)

if(s < datasize + 1){
    logscore <- log(posteriorx(matrix(y[s,], nrow=1),fitdp))
} else {
    othery <- t(as.matrix(read.csv(paste0('weights_', groups[1], '_40cons','.csv'),header=FALSE,sep=',')))
    othery <- othery[,rdims]
    othery <- logit2(othery)
    logscore <- log(posteriorx(othery ,fitdp))
}

write.table(logscore, file= paste0(basename,'_logscore_', group,'_', s, '.csv'), row.names=FALSE, col.names=FALSE, sep=',')

pdff(paste0(basename,'_diagn_',group,'_',s))
DiagnosticPlots(fitdp)
dev.off()


#print('...saved')
#}

#sbatch postdpS.slurm posteriors_S.R
## Idun:
## > dp <- DirichletProcessMvnormal(y,
## +                                 #alphaPriors=c(0.1, 0.1),
## +                                 g0Priors=list(mu0=rep(0,dims), kappa0=0.01, Lambda=diag(dims)/2, nu=dims+1), numInitialClusters = datasize)
## > fitdp <- dp
## > system.time(for(i in 1:1){fitdp <- Fit(fitdp,10000)})
##   |--------------------------------------------------| 100%
##      user    system   elapsed
## 15168.709    47.122   781.522
## > system.time(for(i in 1:1){fitdp <- Fit(fitdp,1000)})
##   |--------------------------------------------------| 100%
##     user   system  elapsed
## 1042.095    3.345   53.594
##
## > system.time(for(i in 1:1){fitdp <- Fit(fitdp,10000)})
##   |--------------------------------------------------| 100%
##      user    system   elapsed
## 10789.939    35.819   558.570
## >


## + system.time(for(i in 1:1){fitdptest <- Fit(fitdptest,5000)})
## + 
##   |--------------------------------------------------| 100%
##    user  system elapsed 
##  317.40    0.15  317.67 

## + system.time(for(i in 1:1){fitdptest <- Fit(fitdptest,1000)})
## + 
##   |--------------------------------------------------| 100%
##    user  system elapsed 
##   52.40    0.15   52.57 


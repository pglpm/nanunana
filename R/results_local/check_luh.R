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

#load('logit-unif_health_all_rdmh.RData')

#tsamples <- logit(samples[,1:d])

for(i in 1:d){
edg[,i] <- seq(quantile(tsamples[,i],0.05), quantile(tsamples[,i],0.95), length.out=nbins+1)
}

pdf(paste0('check_',filename,'_dens.pdf'))
for(i in 1:d){
    posterior <- tsamples[,i]
    dat <- data.frame(x=posterior)
   # data.f <- data.frame(x=xx[i,], y=yy[i,])
    print(
        ggplot() + geom_histogram(data=dat,aes(x=x,y=..density..),
                                  breaks=edg[,i]) #+ scale_y_log10()
        +labs(x=tparm.names[i],title=paste0('mean =',signif(mean(posterior),3),', sd = ',signif(sd(posterior),3)))
        + geom_vline(xintercept=mean(data[,i]),colour=myyellow)
        + geom_vline(xintercept=mean(posterior),colour=mygreen)
        + geom_vline(xintercept=median(posterior),colour=mygreen,linetype='dashed')
        #+ geom_vline(xintercept=datam[,i],colour=myred)
        + geom_point(aes(x=data[,i], y=0),size=2,colour=myred,shape=16)
        #+xlim(-3*sd(posterior),3*sd(posterior))
    )
}
for(j in 1:(d-1)){
    for(i in (j+1):d){
        posterior <- tsamples[,j]
        posterior2 <- tsamples[,i]
        dat <- data.frame(x=posterior,y=posterior2)
        #gmedian <- Gmedian(cbind(posterior,posterior2))
        datapoints <- data.frame(x=data[,j], y=data[,i])
        print(
            ggplot() + stat_bin2d(data=dat,aes(x=x,y=y,fill=..density..),
                                  breaks=list(x=edg[,j], y=edg[,i])
                                  ) +
            scale_fill_gradientn(colours=brewer.pal(n=9,name='Blues')) +
            labs(x=tparm.names[j],y=tparm.names[i],title=ptitle) 
            +geom_point(data=datapoints ,aes(x,y),colour=myred,size=4,shape=16)
            + geom_point(aes(x=mean(posterior), y=mean(posterior2)),size=3,colour=mygreen,shape=17)
            + geom_point(aes(x=mean(data[,j]), y=mean(data[,i])),size=3,colour=myyellow,shape=15)
            #+ geom_point(aes(x=gmedian[1], y=gmedian[2]),size=3,colour=mygreen,shape=15)
            #+ coord_fixed()
            #+xlim(-3*sd(posterior),3*sd(posterior))
            #+ylim(-3*sd(posterior2),3*sd(posterior2))
        )
        }
}
dev.off()    

library('fBasics')
library('LaplacesDemon')
library('ggplot2')
library('RColorBrewer')
library('mvtnorm')
##library('Matrix')
mypurpleblue <- '#4477AA'
myblue <- '#66CCEE'
mygreen <- '#228833'
myyellow <- '#CCBB44'
myred <- '#EE6677'
myredpurple <- '#AA3377'
mygrey <- '#BBBBBB'
palette(c(myblue, myred, mygreen, myyellow, myredpurple, mypurpleblue, mygrey, 'black'))
dev.off()


jlog <- function(x){1/x}
jlogit <- function(x){1/(x*(1-x))}

mvn <- function(x,mu,sigma){
    exp(-c((x - mu) %*% solve(sigma) %*% (x - mu))/2)/sqrt(det(2*pi*sigma))
}


likelihood <- function(x,parameters){
    l <- length(x)
    kappa <- parameters[[2]]
    nu <- parameters[[1]] - l + 1
    dmvt(x, delta=parameters[[3]],
         sigma=parameters[[4]]*(kappa+1)/kappa/nu, df=nu, log=F)
}
loglikelihood <- function(x,parameters){
    l <- length(x)
    kappa <- parameters[[2]]
    nu <- parameters[[1]] - l + 1
    dmvt(x, delta=parameters[[3]],
         sigma=parameters[[4]]*(kappa+1)/kappa/nu, df=nu, log=T)
}

updateparameters <- function(data,parameters){
    N <- dim(data)[2]
    means <- colMeans(t(data))
    scatm <- cov(t(data))*(N-1)
    nu <- parameters[[1]]
    kappa <- parameters[[2]]
    newkappa <- kappa + N
    mu <- parameters[[3]]
    sigma <- parameters[[4]]
    list(nu + N, newkappa,
    (kappa*mu + N * means)/newkappa,
    sigma + scatm + kappa * N * ((means-mu) %*% t(means-mu))/newkappa )
} # checked against same alg in mathematica
updateparameters1 <- function(datum,parameters){
    nu <- parameters[[1]]
    kappa <- parameters[[2]]
    newkappa <- kappa + 1
    mu <- parameters[[3]]
    sigma <- parameters[[4]]
    list(nu + 1, newkappa,
    (kappa*mu + datum)/newkappa,
    sigma + kappa * ((datum-mu) %*% t(datum-mu))/newkappa )
} # checked against same alg in mathematica

dd <- 25
priorparameters0d <- function(d){list(d+1, 1, rep(0,d), (2+d)*diag(rep(1,d)))}
priorparameters0 <- priorparameters0d(dd)
priorparametersf <- function(d,kappa,mu){list(d+1, kappa, rep(0,d), mu*(2+d)*diag(rep(1,d)))}
cupprior <- list(dd + 1, 0.1, rep(0,dd), (2*dd+3)*diag(rep(1,dd)))


## mvt <- function(x,nu,mu,sigma){
##     l <- length(mu)
##    Gamma[(nu + l)/2]/Gamma[nu/2]/
##      Sqrt@Det[
##        nu*Pi*sigma]/(1 + (x - mu).Inverse[nu*sigma].(x - mu))^((nu + 
##          l)/2)];
## }

graphprop <- c("graphweights_1", "graphweights_2", "graphweights_3", "graphweights_4", "shortest path_1", "shortest path_2", "shortest path_3", "shortest path_4", "cluster coeff_1", "cluster coeff_2", "cluster coeff_3", "cluster coeff_4", "weighted degree_1", "weighted degree_2", "weighted degree_3", "weighted degree_4", "norm. weighted degree_1", "norm. weighted degree_2", "norm. weighted degree_3", "norm. weighted degree_4", "closeness centrality1", "closeness centrality2", "closeness centrality3", "closeness centrality4", "node size")
transf <- c(logit, logit, logit, logit, log, log, log, log, logit, logit, logit,logit, log, log, log, log, log, log, log, log, logit, logit, logit, logit, log)
transfd <- function(data){mapply(function(f,x) f(x),transf,data)}
scales <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 100, 100, 100, 100, 1, 1, 1, 1, 1, 1, 1, 1, 100)

logevidence <- function(nsamples,categories,quantities,datafiles,prior,pretestprob,seed=666){
    n <- length(categories) # num. health categories
    d <- dim(read.matrix(datafiles[1]))[1] # num. graph quantities
    if(length(quantities)==1){quantities <- sprintf("q[%d]",1:d)}
    datao <- list() # original data
    datas <- list() # scaled data
    datat <- list() # logit/log-transformed scaled data
    N <- rep(NA,n) # num. patients per health category
    for(i in 1:n){
        categ <- categories[i]
        ## read data
        datao[[categ]] <- read.matrix(datafiles[i])#[1:d,1:(N[i])]
        N[i] <- dim(datao[[categ]])[2]
        ## replace central moments with raw moments
        tempdata <- datao[[categ]]
        datanew <- datao[[categ]]
        ## transform 2nd moments to raw ones
        for(j in seq(1,21,4)){
            datanew[j + 1,] <- sqrt(tempdata[j,]^2 + tempdata[j + 1,])
            ## transform 3nd moments to raw ones
            datanew[j + 2,] <- (tempdata[j,]^3 + 3*tempdata[j,]*tempdata[j + 1,] + tempdata[j + 2,])^(1/3)
            ## transform 4nd moments to raw ones
            tempdata[j + 3,] <- (tempdata[j,]^4 + 6*tempdata[j,]^2*tempdata[j + 1,] + 4*tempdata[j,]*tempdata[j + 2,] + tempdata[j + 3,])^(1/4)
        }
        datao[[categ]] <- datanew
        
        rownames(datao[[categ]]) <- quantities
        colnames(datao[[categ]]) <- sprintf("id[%d]",1:(N[i]))
        datas[[categ]] <- datao[[categ]]/scales # scaled data
        ## create logit/log-transformed scaled data
        datat[[categ]] <- datas[[categ]]
        for(j in 1:(N[i])){
           datat[[categ]][,j] <- mapply(function(f,x) f(x),transf,datas[[categ]][,j])
        }
    }
    NN <- sum(N) # total number of patients
    lmlsamples1 <- matrix(NA,nsamples,NN-1) # log-evidence from graph quantities
    lmlsamples2 <- lmlsamples1 # log-evidence from health condition
    probsamples1 <- matrix(NA,nsamples,NN-1) # graph quantity density for next patient
    logprobsamples1 <- probsamples1 # log of above
    probsamples2 <- probsamples1 # health probability for next patient

    for(k in 1:nsamples){set.seed(seed+k) # if we want the same patient sequence
        train <- sample(NN) # randomize order of patients
        lml1 <- 0
        lml2 <- 0
        ## set parameters to initial prior
        newprior <- list()
        for(i in 1:n){newprior[[i]] <- prior}
        
        for(i in 1:(NN-1)){ # take in randomized training data in succession
            ## loop to find out health category of current patient
            patient <- train[i] # current patient
            categ <- 1
            while(patient > N[categ]){
                patient <- patient - N[categ]
                categ <- categ + 1}
            datum <- datat[[categ]][,patient] # graph data for current patient

            ## log-evidence based on graph quantities
            prob1 <- likelihood(datum,newprior[[categ]])
            logprob1 <- loglikelihood(datum,newprior[[categ]])
            lml1 <- lml1 + logprob1

            ## log-evidence based on health condition
            ## normalization factor for probability of health condition
            sumprobs <- prob1 * pretestprob[categ]
            for(j in (1:n)[-categ]){
                sumprobs <- sumprobs +
                    likelihood(datum,newprior[[j]])*pretestprob[j]}
            lml2 <- lml2 + logprob1 + log(pretestprob[categ]) - log(sumprobs)

            
            ## prediction about next patient (some variables are reassigned here)
            ## update parameters with current patient
            newprior[[categ]] <- updateparameters1(datum,newprior[[categ]])
            ## loop to find out health category of next patient
            patient <- train[i+1] # next patient
            categ <- 1
            while(patient > N[categ]){
                patient <- patient - N[categ]
                categ <- categ + 1}
            datum <- datat[[categ]][,patient] # graph data for next patient
            prob1 <- likelihood(datum,newprior[[categ]])
            logprob1 <- loglikelihood(datum,newprior[[categ]])
            
            prob2 <- prob1 * pretestprob[categ]
            ## normalization factor for probability of health condition
            sumprobs <- prob2
            for(j in (1:n)[-categ]){
                sumprobs <- sumprobs +
                    likelihood(datum,newprior[[j]])*pretestprob[j]}
            prob2 <- prob2/sumprobs

            lmlsamples1[k,i] <- lml1
            lmlsamples2[k,i] <- lml2
            probsamples1[k,i] <- prob1
            logprobsamples1[k,i] <- logprob1
            probsamples2[k,i] <- prob2
        }
    }
    return(list(lmlsamples1=lmlsamples1,
                lmlsamples2=lmlsamples2,
                probsamples1=probsamples1,
                logprobsamples1=logprobsamples1,
                probsamples2=probsamples2))
}

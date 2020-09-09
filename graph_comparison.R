## comparison of predictions of health conditions by different graph analyses

## file 'logitnormal_evidence.R' needs to be in same directory
source('logitnormal_evidence.R')

## For each graph analysis we need a list with the following:
##
## 1-3. Files with graph values for healthy, Alzheimer, MCI patients. One column for each graph property, one row for each patient
##
## 4. Vector of the form c(logit, logit, log, log, ...), with one entry for each graph property: logit if the graph property is in [0,1], log if the graph property is in [0,inf). The order is the same as the rows of the data files for the patients.
##
## 5. Vector of the form c(jlogit, jlogit, jlog, jlog, ...),  with one entry for each graph property: jlogit if the graph property is in [0,1], jlog if the graph property is in [0,inf). These are the Jacobians for the changes of variable.
##
## 6. Vector of the form c(1, 1, 1, 100, ...), with one entry for each graph property. The values are the rescalings of the respective graph properties. [0,1] quantities should all have scale 1. [0,inf) quantities may have scale 1, or 10, or 100, or 0.1, etc. The goal is that the averages of their log quantities should be as close to 1 as possible.
##
## These are the functions and scales for the "richclub" quantities:
##
## "graphweights_1"		logit		1		
## "graphweights_2"		logit		1		
## "graphweights_3"		logit		1		
## "graphweights_4"		logit		1		
## "shortest path_1"		log		1		
## "shortest path_2"		log		1		
## "shortest path_3"		log		1		
## "shortest path_4"		log		1		
## "cluster coeff_1"		logit		1		
## "cluster coeff_2"		logit		1		
## "cluster coeff_3"		logit		1		
## "cluster coeff_4"		logit		1		
## "weighted degree_1"		log		100
## "weighted degree_2"		log		100
## "weighted degree_3"		log		100
## "weighted degree_4"		log		100
## "norm. weighted degree_1"	log		1		
## "norm. weighted degree_2"	log		1		
## "norm. weighted degree_3"	log		1		
## "norm. weighted degree_4"	log		1		
## "closeness centrality1"	logit		1		
## "closeness centrality2"	logit		1		
## "closeness centrality3"	logit		1		
## "closeness centrality4"	logit		1		
## "node size"			log		100



## Example

healthconditions <- c('H', 'AD', 'MCI')

graphlist <- list()

## richclub
graphlist[[1]] <- list(
    c('Control0_8', 'Alzheimer0_8', 'MCI0_8'),
    c(log, log, log, log, log, log, log, log, logit, logit, logit, logit, log, log, log, log, logit, logit, logit, logit, log, log),
    c(jlog, jlog, jlog, jlog, jlog, jlog, jlog, jlog, jlogit, jlogit, jlogit, jlogit, jlog, jlog, jlog, jlog, jlogit, jlogit, jlogit, jlogit, jlog, jlog),
    c(1.0, 1.0, 1.0, 1.0, 0.01, 0.01, 0.01, 0.01, 1, 1, 1, 1, 0.1, 0.1, 0.1, 0.1, 1, 1, 1, 1, 1.0, 100.0)
)



numberofanalyses <- length(graphlist)
numberofsamples <- 560 ## 100 gives rougher results but is much faster
prior <- priorparameters0 ## the one that's worked best so far
pretestprob <- c(1,1,1)/3

## this calculates the log-evidence and "hits" for each graph analysis in turn
logevidencelist <- list()
hitslist <- list()
for(i in 1:numberofanalyses){
print(paste0('Calculating graph setup ',i))
    filenames <- graphlist[[i]][[1]]
    transformations <- graphlist[[i]][[2]]
    scales <- graphlist[[i]][[4]]
    jacobians <- graphlist[[i]][[3]]
    
    results <- logevidence(numberofsamples,healthconditions,c(1),
                           filenames,
                           transformations,
                           scales,
                           jacobians,
                           prior,
                           pretestprob)
    logevidencelist[[i]] <- prob2meanlogev(results)
    hitslist[[i]] <- prob2meanhits(results)
}

## transforms the resulting lists to matrices
logevidencedata <- matrix(unlist(logevidencelist),
                              nrow=numberofanalyses,byrow=T)
hitsdata <- matrix(unlist(hitslist),
                       nrow=numberofanalyses,byrow=T)

## Then we can plot the evolution of log-evidence with training for all graph analyses

print('Plotting results...')

logevpdfname <- 'results_logevidence13.pdf'
hitspdfname <- 'results_hits13.pdf'

meanlogevplot(2, logevidencedata,logevpdfname,yrange=NULL) ## y range can be set to something else if needed

meanhitsplot(hitsdata,hitspdfname,yrange=NULL) ## y range can be set to something else if needed
write.csv(logevidencedata, file='le_data')
write.csv(hitsdata, file='hitsdata' )
print(paste0('Final log-evidences: ',logevidencedata[,dim(logevidencedata)[2]]))
print('')
print(paste0('Final hits (leave-1-out cross-val): ',hitsdata[,dim(hitsdata)[2]]))


############################################################
## Implement the Reed-McCann Method (independent case)
# Coded by Taotao Zhou, Ph.D., 
# College of Safety and Ocean Engineering, 
# China University of Petroleum-Beijing, Beijing, China
# Please cite:
# (1) T. Zhou, M. Modarres and E.L. Droguett. “Issues in Dependency Modeling 
# in Multi- Unit Seismic PRA.” (PSA2017),September 24-28, 2017, Pittsburgh, PA.
# (2) T. Zhou, M. Modarres and E.L. Droguett. “An Improved Multi-Unit Nuclear 
# Plant Seismic Probabilistic Risk Assessment Approach.” Reliability Engineering 
# & System Safety, 171 (2018), 34-47.
# (3) T. Zhou, L. Zhang, J. Hu, M. Modarres and E.L. Droguett. “A Critical 
# Review and Benchmark Study of Dependency Modeling for Seismic Probabilistic 
# Risk Assessment in the Nuclear Power Industry.” Reliability Engineering & 
# System Safety (2024), 110009.    
############################################################
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!

#------------------------ Initialize the Working Path ------------------------#
mydirectory <- 'C:\\Users\\taota\\Desktop\\PSAM_seismic\\Reed-McCann'
setwd(mydirectory)   # change to mydirectory
options(digits = 12)

############################################################
## Latin Hypercube Sampling
pr = c(0.02,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.98)
stdvar = qnorm(pr)
num.sample = length(pr)-1 # number of samples

## Input Fragility
Am = c(0.9,1.0,1.1)
BetaR = c(0.6,0.5,0.6)
BetaU = c(0.6,0.5,0.4)
BetaC = sqrt(BetaR^2+BetaU^2)

## Common Beta Values
BetaRS = c(0.3,0.45,0)
BetaUS = c(0.4,0.35,0)
BetaCS = sqrt(BetaRS^2+BetaUS^2)

## Reduced Beta Values
BetaRP = c(sqrt(BetaR[1]^2-BetaRS[1]^2-BetaRS[2]^2),
           sqrt(BetaR[2]^2-BetaRS[1]^2-BetaRS[3]^2),
           sqrt(BetaR[3]^2-BetaRS[2]^2-BetaRS[3]^2))
BetaUP = c(sqrt(BetaU[1]^2-BetaUS[1]^2-BetaUS[2]^2),
           sqrt(BetaU[2]^2-BetaUS[1]^2-BetaUS[3]^2),
           sqrt(BetaU[3]^2-BetaUS[2]^2-BetaUS[3]^2))


############################################################
## Stage 1
############################################################
## Generate LN Ordinates
varblps = matrix(0,nrow=length(pr),ncol=6)
for (i in 1:3){
  varblps[,i] = Am[i]*exp(BetaUP[i]*stdvar)
  
}
for (i in 4:6){
  varblps[,i] = 1*exp(BetaUS[i-3]*stdvar)
}

## Ten random samples
sampleABCDE = matrix(0,nrow=num.sample,ncol=6)
for (i in 1:6){
  for (j in 1:num.sample){
    sampleABCDE[j,i] = runif(1,varblps[j,i],varblps[j+1,i])
  }
}

## Shuffled samples
ranperm <- function(X, N){
  order(runif(N))
} 
locABC = matrix(nrow = num.sample, ncol = 6)
locABC = apply(locABC, 2, ranperm, N = num.sample)

shuffledABCDE = matrix(0,nrow=num.sample,ncol=6)
for (i in 1:ncol(sampleABCDE))
{
  index = locABC[,i]
  shuffledABCDE[index,i] = sampleABCDE[,i]
}

## Combine values
XABCc = matrix(0,nrow=num.sample,ncol=3)
XABCc[,1] = shuffledABCDE[,1]*shuffledABCDE[,4]*shuffledABCDE[,5]
XABCc[,2] = shuffledABCDE[,2]*shuffledABCDE[,4]*shuffledABCDE[,6]
XABCc[,3] = shuffledABCDE[,3]*shuffledABCDE[,5]*shuffledABCDE[,6] 

############################################################
## Stage 2
############################################################
## Enter discrete acceleration values
# agn = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)
agn = c(0.25,0.45,0.65,0.85,1.00,1.10,1.20,1.30,1.40,1.50)

############################################################
## Multiple Integration Approach
############################################################
## Define variables used in sampling and triple integral

outInterMI = matrix(0,nrow=length(agn),ncol=num.sample)
outUnionMI = matrix(0,nrow=length(agn),ncol=num.sample)


for (ipga in 1:length(agn))
{
  a = agn[ipga]
  
  for (iobs in 1:num.sample)
  {
    ## Intersection
    fInter <- function(X) {
      1/(BetaRS[1]*BetaRS[2])*
        pnorm(log(a/X[1]/X[2]/XABCc[iobs,1])/BetaRP[1])*
        pnorm(log(a/X[1]/XABCc[iobs,2])/BetaRP[2])*
        pnorm(log(a/X[2]/XABCc[iobs,3])/BetaRP[3])*
        dnorm(log(X[1])/BetaRS[1])*
        dnorm(log(X[2])/BetaRS[2])/
        (X[1]*X[2])
    } # "x" is vector
    
    ## Union
    fUnion <- function(X) {
      1/(BetaRS[1]*BetaRS[2])*
        (1-(1-pnorm(log(a/X[1]/X[2]/XABCc[iobs,1])/BetaRP[1]))*
           (1-pnorm(log(a/X[1]/XABCc[iobs,2])/BetaRP[2]))*
           (1-pnorm(log(a/X[2]/XABCc[iobs,3])/BetaRP[3])))*
        dnorm(log(X[1])/BetaRS[1])*
        dnorm(log(X[2])/BetaRS[2])/
        (X[1]*X[2])
    } # "x" is vector
    
    library(cubature) # load the package "cubature"
    outInter = adaptIntegrate(fInter, lowerLimit = c(0, 0), upperLimit = c(2.1, 2.1))
    outInterMI[ipga,iobs]=outInter$integral
    outUnion = adaptIntegrate(fUnion, lowerLimit = c(0, 0), upperLimit = c(2.1, 2.1))
    outUnionMI[ipga,iobs]=outUnion$integral
    
  }
  
}

########## Results Summary
outputInterMI = cbind(agn,outInterMI,rowMeans(outInterMI)) #parallel system
outputUnionMI = cbind(agn,outUnionMI,rowMeans(outUnionMI)) #serial system

########## Calculate failure probability
prob_parallel = rowMeans(outInterMI) #parallel system
prob_serial = rowMeans(outputUnionMI) #serial system

########## output results
write.table(prob_parallel,"Reed_parallel.csv", row.names = FALSE, col.names = FALSE, sep=',')

write.table(prob_serial,"Reed_serial.csv", row.names = FALSE, col.names = FALSE, sep=',')

########## Visualization
par(mfrow=c(1,2))
## Plot for intersection
plot(prob_parallel,type = "l", col = 'red',lty = 1:5, lwd = 1,xlim = NULL, ylim = NULL,
     xlab = 'Acceleration(g)', ylab = 'Conditional Failure Probability')
points(prob_parallel,lwd = 3,col = 'red',pch=21)
title(main = 'Failure Probability(Parallel)', sub = NULL, cex.main = 0.8)

plot(prob_serial,type = "l", col = 'red',lty = 1:5, lwd = 1,xlim = NULL, ylim = NULL,
     xlab = 'Acceleration(g)', ylab = 'Conditional Failure Probability')
points(prob_serial,lwd = 3,col = 'red',pch=21)
title(main = 'Failure Probability(Serial)', sub = NULL, cex.main = 0.8)



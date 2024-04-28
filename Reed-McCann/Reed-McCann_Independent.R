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
## Ten random samples
varblps.Ind = matrix(0,nrow=length(pr),ncol=3)
for (i in 1:3){
  varblps.Ind[,i] = Am[i]*exp(BetaU[i]*stdvar)
  
}

## Ten random samples
sampleABCDE.Ind = matrix(0,nrow=num.sample,ncol=3)
for (i in 1:3){
  for (j in 1:num.sample){
    sampleABCDE.Ind[j,i] = runif(1,varblps.Ind[j,i],varblps.Ind[j+1,i])
  }
}

## Shuffled samples
ranperm <- function(X, N){
  order(runif(N))
} 
locABC.Ind = matrix(nrow = num.sample, ncol = 3)
locABC.Ind = apply(locABC.Ind, 2, ranperm, N = num.sample)

shuffledABCDE.ind = matrix(0,nrow=num.sample,ncol=3)
for (i in 1:ncol(sampleABCDE.Ind))
{
  index = locABC.Ind[,i]
  shuffledABCDE.ind[index,i] = sampleABCDE.Ind[,i]
}


############################################################
## Stage 2
############################################################
## Enter discrete acceleration values
# agn = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)
agn = c(0.25,0.45,0.65,0.85,1.00,1.10,1.20,1.30,1.40,1.50)

#################################################
## Multiple Integration Approach
############################################################
## Define variables used in sampling and triple integral

# Sample Set of Interest: no.1

outInterIndMI = matrix(0,nrow=length(agn),ncol=num.sample)
outUnionIndMI = matrix(0,nrow=length(agn),ncol=num.sample)

for (ipga in 1:length(agn))
{
  a = agn[ipga]
  
  for (iobs in 1:num.sample)
  {
    ###################
    # Independent Case
    ## Intersection
    fInterInd <- exp(pnorm(log(a/shuffledABCDE.ind[iobs,1])/BetaR[1],log.p = TRUE))*
      exp(pnorm(log(a/shuffledABCDE.ind[iobs,2])/BetaR[2],log.p = TRUE))*
      exp(pnorm(log(a/shuffledABCDE.ind[iobs,3])/BetaR[3],log.p = TRUE))
    
    ## Union
    fUnionInd <- 1-pnorm(log(a/shuffledABCDE.ind[iobs,1])/BetaR[1],lower.tail = FALSE)*
      pnorm(log(a/shuffledABCDE.ind[iobs,2])/BetaR[2],lower.tail = FALSE)*
      pnorm(log(a/shuffledABCDE.ind[iobs,3])/BetaR[3],lower.tail = FALSE)
    
    outInterIndMI[ipga,iobs] = fInterInd
    outUnionIndMI[ipga,iobs] = fUnionInd
    
  }
  
}

########## Results Summary
outputInterIndMI = cbind(agn,outInterIndMI,rowMeans(outInterIndMI))
outputUnionIndMI = cbind(agn,outUnionIndMI,rowMeans(outUnionIndMI))

########## Calculate failure probability
prob_parallel_ind = rowMeans(outputInterIndMI) #parallel system
prob_serial_ind = rowMeans(outputUnionIndMI) #serial system


########## output results
write.table(prob_parallel_ind,"Reed_parallel_ind.csv", row.names = FALSE, col.names = FALSE, sep=',')

write.table(prob_serial_ind,"Reed_serial_ind.csv", row.names = FALSE, col.names = FALSE, sep=',')

########## Visualization
par(mfrow=c(1,2))
## Plot for intersection
plot(prob_parallel_ind,type = "l", col = 'red',lty = 1:5, lwd = 1,xlim = NULL, ylim = NULL,
     xlab = 'Acceleration(g)', ylab = 'Conditional Failure Probability')
points(prob_parallel_ind,lwd = 3,col = 'red',pch=21)
title(main = 'Failure Probability(Parallel_IND)', sub = NULL, cex.main = 0.8)

plot(prob_serial_ind,type = "l", col = 'red',lty = 1:5, lwd = 1,xlim = NULL, ylim = NULL,
     xlab = 'Acceleration(g)', ylab = 'Conditional Failure Probability')
points(prob_serial_ind,lwd = 3,col = 'red',pch=21)
title(main = 'Failure Probability(Serial_IND)', sub = NULL, cex.main = 0.8)

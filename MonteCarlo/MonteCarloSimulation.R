############################################################
## Implement the Monte Carlo Simulation Method (independent and dependent cases)
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
mydirectory <- 'C:\\Users\\taota\\Desktop\\PSAM_seismic\\MonteCarlo'
setwd(mydirectory)   # change to mydirectory
options(digits = 12)

library(msm)
library(copula)

#------------------------ Discretize Ground Motion Intensity ------------------------#
# Conditional on ground motion
smotion <- c(0.05,0.25,0.45,0.65,0.85,1,1.1,1.2,1.3,1.4,1.5)
disGroundMotion = replicate((length(smotion)), list()) 
for (imotion in 1:(length(smotion)-1)) {
  disGroundMotion[[imotion]] <- list(smin = smotion[imotion], smax=smotion[imotion+1] )
}
num.pgaInterval = length(smotion)-1

#------------------------ Set up importance sampling ------------------------#
# Target distribution
p <- function(t) (1*10^-5*t^-1.61) # with normalization constanct as 6*10^-6/0.00618 = 9.70874e-4
## Proposal distribution - Truncated Normal
q <- function(t,mean,sd,tmin,tmax) (dtnorm(t,mean,sd,lower=tmin, upper=tmax))
## Indicator Function: Binary State for Components
f <- function(cf,t){
  cs <- as.numeric(cf <= t)
  return(cs) # component states failure or survive
}

#------------------------ System INPUT - Seismic Failure ------------------------#
# component_1
Capacity_1 = c(0.9, sqrt(0.6^2+0.6^2))
# component_2
Capacity_2 = c(1.0, sqrt(0.5^2+0.5^2))
# component_3
Capacity_3 = c(1.1, sqrt(0.4^2+0.6^2))

#------------------------ Specify Paramters for Monte Carlo ------------------------#
# num.iterations = 10  # no. of monte carol runs to execute
samp.size <- 1e6 # no. of initial samples to draw from sampling distribution
# Specify the OUTPUT
outputFinal = matrix(0,nrow=num.pgaInterval,ncol=4) 

#------------------------ Execute Monte Carlo ------------------------#
for (ipga in 1:num.pgaInterval)
{
  tmin = disGroundMotion[[ipga]]$smin
  tmax= disGroundMotion[[ipga]]$smax

  ########################################################################
  ## Execute Importance Sampling
  mean.proposal = tmax
  sd.proposal = (tmax-tmin)/1
  # Generate ground motion intensity from proposal distribution
  X <- rtnorm(samp.size,mean.proposal,sd.proposal,lower=tmin, upper=tmax)
  # Check the samples from proposal distribution
  hist(X)
  # Determine the weight vector for samples
  W <- p(X)/q(X,mean.proposal,sd.proposal,tmin,tmax) # calculate importance weights for samples
  
  ########################################################################
  ## Generate independent capacity
  # component_1
  indCap_1 = rlnorm (samp.size, meanlog = log(Capacity_1[1]),
                     sdlog = Capacity_1[2])
  # component_2
  indCap_2 = rlnorm (samp.size, meanlog = log(Capacity_2[1]),
                     sdlog = Capacity_2[2])
  # component_3
  indCap_3 = rlnorm (samp.size, meanlog = log(Capacity_3[1]),
                     sdlog = Capacity_3[2])
  
  ## Determine component states - Independent Case
  state_1_IND = as.numeric(indCap_1<=X)
  state_2_IND = as.numeric(indCap_2<=X)
  state_3_IND = as.numeric(indCap_3<=X)
  
  ## determine system state for independent(parallel system)
  state_Ind_parallel = state_1_IND*state_2_IND*state_3_IND
  cdf_Ind_parallel= sum(W*as.numeric(state_Ind_parallel>=1))/sum(W)
  
  ## determine system state for independent(serial system)
  state_Ind_serial = state_1_IND+state_2_IND+state_3_IND
  cdf_Ind_serial= sum(W*as.numeric(state_Ind_serial>=1))/sum(W)

  ########################################################################
  ## Generate dependent capacity and replace the capacities of Correlated Component
  # Construct Normal copula with Normal margins
  myCop <- normalCopula(c(0.42, 0.53, 0), dim = 3, dispstr = "un")
  myMvd <- mvdc(copula=myCop, margins=rep("lnorm",3), 
                paramMargins= list(list(meanlog=log(Capacity_1[1]), sdlog=Capacity_1[2]), 
                                   list(meanlog=log(Capacity_2[1]), sdlog=Capacity_2[2]), 
                                   list(meanlog=log(Capacity_3[1]), sdlog=Capacity_3[2])))
  # Generate correlated random variables from copula
  output <- rMvdc(samp.size,myMvd)
  # component_1_dep
  depCap_1 = matrix(output[,1])
  # component_2_dep
  depCap_2 = matrix(output[,2])
  # component_3_dep
  depCap_3 = matrix(output[,3])
  
  ## Determine component states - Dependent Case
  state_1_DEP = as.numeric(depCap_1<=X)
  state_2_DEP = as.numeric(depCap_2<=X)
  state_3_DEP = as.numeric(depCap_3<=X)
  SeismicCompState = rowSums(cbind(state_1_DEP,state_2_DEP,state_3_DEP))
  
  ## determine system state for dependent (parallel system)
  state_Dep_parallel = state_1_DEP*state_2_DEP*state_3_DEP
  cdf_Dep_parallel= sum(W*as.numeric(state_Dep_parallel>=1))/sum(W)
  
  ## determine system state for dependent (serial system)
  state_Dep_serial = state_1_DEP+state_2_DEP+state_3_DEP
  cdf_Dep_serial = sum(W*as.numeric(state_Dep_serial>=1))/sum(W)
  
  ########################################################################
  ## output
  #ind_serial,ind_parallel,dep_serial,dep_parallel
  outputFinal[ipga,] = c(cdf_Ind_serial, cdf_Ind_parallel, cdf_Dep_serial, cdf_Dep_parallel)
  
}


########## output results
write.table(outputFinal[,1],"MCS_ind_serial.csv", row.names = FALSE, col.names = FALSE, sep=',')

write.table(outputFinal[,2],"MCS_ind_parallel.csv", row.names = FALSE, col.names = FALSE, sep=',')

write.table(outputFinal[,3],"MCS_serial.csv", row.names = FALSE, col.names = FALSE, sep=',')

write.table(outputFinal[,4],"MCS_parallel.csv", row.names = FALSE, col.names = FALSE, sep=',')

########## Visualization-ind
par(mfrow=c(1,2))
## Plot for intersection
plot(outputFinal[,2],type = "l", col = 'red',lty = 1:5, lwd = 1,xlim = NULL, ylim = NULL,
     xlab = 'Acceleration(g)', ylab = 'Conditional Failure Probability')
points(outputFinal[,2],lwd = 3,col = 'red',pch=21)
title(main = 'Failure Probability(Parallel_IND)', sub = NULL, cex.main = 0.8)

plot(outputFinal[,1],type = "l", col = 'red',lty = 1:5, lwd = 1,xlim = NULL, ylim = NULL,
     xlab = 'Acceleration(g)', ylab = 'Conditional Failure Probability')
points(outputFinal[,1],lwd = 3,col = 'red',pch=21)
title(main = 'Failure Probability(Serial_IND)', sub = NULL, cex.main = 0.8)

########## Visualization-dep
par(mfrow=c(1,2))
## Plot for intersection
plot(outputFinal[,4],type = "l", col = 'red',lty = 1:5, lwd = 1,xlim = NULL, ylim = NULL,
     xlab = 'Acceleration(g)', ylab = 'Conditional Failure Probability')
points(outputFinal[,4],lwd = 3,col = 'red',pch=21)
title(main = 'Failure Probability(Parallel)', sub = NULL, cex.main = 0.8)

plot(outputFinal[,3],type = "l", col = 'red',lty = 1:5, lwd = 1,xlim = NULL, ylim = NULL,
     xlab = 'Acceleration(g)', ylab = 'Conditional Failure Probability')
points(outputFinal[,3],lwd = 3,col = 'red',pch=21)
title(main = 'Failure Probability(Serial)', sub = NULL, cex.main = 0.8)







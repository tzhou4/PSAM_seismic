############################################################
## Implement the COREX Method (dependent case)
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
mydirectory <- 'C:\\Users\\taota\\Desktop\\Software\\COREX'

setwd(mydirectory)   # change to mydirectory
options(digits = 12)
ptm <- proc.time()

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
# q <- function(t,mean,sd) (dtnorm(t,mean,sd,lower=tmin, upper=tmax))
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
outputFinal = matrix(0,nrow=num.pgaInterval,ncol=9) 
#ind,dep,beta

#------------------------ Execute Monte Carlo ------------------------#
for (ipga in 1:num.pgaInterval)
{
  tmin = disGroundMotion[[ipga]]$smin
  tmax= disGroundMotion[[ipga]]$smax
  #if don't consider uncertainty in discrete interval
  #tref = (tmin+tmax)/2; #tref = sqrt(tmin*tmax); #tref = tmax
  
  ########################################################################
  ## Execute Importance Sampling
  #mean.proposal = sqrt(tmin*tmax)
  mean.proposal = tmax
  sd.proposal = (tmax-tmin)/1
  # Generate ground motion intensity from proposal distribution
  X <- rtnorm(samp.size,mean.proposal,sd.proposal,lower=tmin, upper=tmax)
  # Check the samples from proposal distribution
  hist(X)
  # Determine the weight vector for samples
  # W <- p(X)/q(X,mean.proposal,sd.proposal) # calculate importance weights for samples
  W <- p(X)/q(X,mean.proposal,sd.proposal,tmin,tmax)
  
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
  
  p1 = sum(W*as.numeric(state_1_DEP>=1))/sum(W)
  p2 = sum(W*as.numeric(state_2_DEP>=1))/sum(W)
  p3 = sum(W*as.numeric(state_3_DEP>=1))/sum(W)
  
  s12 = state_1_DEP+state_2_DEP
  p12 = sum(W*as.numeric(s12>=1))/sum(W)
  
  s13 = state_1_DEP+state_3_DEP
  p13 = sum(W*as.numeric(s13>=1))/sum(W)
  
  s23 = state_2_DEP+state_3_DEP
  p23 = sum(W*as.numeric(s23>=1))/sum(W)
  
  s123 = state_1_DEP+state_2_DEP+state_3_DEP
  p123 = sum(W*as.numeric(s123>=1))/sum(W)
  
  yy = c(p1, p2, p3, p12, p13, p23, p123)
  
  ##solve ccf element
  A <- c(1,0,0,1,1,0,1,
         0,1,0,1,0,1,1,
         0,0,1,0,1,1,1,
         1,1,0,1,1,1,1,
         1,0,1,1,1,1,1,
         0,1,1,1,1,1,1,
         1,1,1,1,1,1,1)
  A <- matrix(A, nrow = 7, byrow = TRUE)
  
  B <- log(1-yy)
  
  sol <- solve(A, B)
  
  X <- 1-exp(sol)
  
  Q1 <- X[1]
  Q2 <- X[2]
  Q3 <- X[3]
  Q12 <- X[4]
  Q13 <- X[5]
  Q23 <- X[6]
  Q123 <- X[7]
  
  parallel_S = Q1*Q2*Q3 + Q123 + Q1*Q23 + Q2*Q13 + Q3*Q12 + Q12*Q13 + Q13*Q23 + Q12*Q23 
  serial_S = 1-(1-Q1)*(1-Q2)*(1-Q3)*(1-Q12)*(1-Q13)*(1-Q23)*(1-Q123)
    
  outputFinal[ipga,] = c(X,parallel_S,serial_S)
}


#Total processing time
timeUse = proc.time() - ptm

outputFinal

plot (outputFinal[,8]) #PARALLEL
plot (outputFinal[,9]) #SERIAL

########## output
write.table(outputFinal[,9],"COREX_serial.csv", row.names = FALSE, col.names = FALSE, sep=',')

write.table(outputFinal[,8],"COREX_parallel.csv", row.names = FALSE, col.names = FALSE, sep=',')


# y1=outputFinal[,1]
# y2=outputFinal[,2]
# y3=outputFinal[,3]
# y12=outputFinal[,4]
# y13=outputFinal[,5]
# y23=outputFinal[,6]
# y123=outputFinal[,7]
# 
# 
# A <- c(1,0,0,1,1,0,1,
#       0,1,0,1,0,1,1,
#       0,0,1,0,1,1,1,
#       1,1,0,1,1,1,1,
#       1,0,1,1,1,1,1,
#       0,1,1,1,1,1,1,
#       1,1,1,1,1,1,1)
# A <- matrix(A, nrow = 7, byrow = TRUE)
# 
# B <- log(1-outputFinal[8,])
# 
# 
# sol <- solve(A, B)
# 
# X <- 1-exp(sol)
# 
# Q1 <- X[1]
# Q2 <- X[2]
# Q3 <- X[3]
# Q12 <- X[4]
# Q13 <- X[5]
# Q23 <- X[6]
# Q123 <- X[7]
# 
# parallel_S = Q1*Q2*Q3 + Q123 + Q1*Q23 + Q2*Q13 + Q3*Q12 + Q12*Q13 + Q13*Q23 + Q12*Q23 

# x1 = log(1-Q1)
# x2 = log(1-Q2)
# x3 = log(1-Q3)
# x12 = log(1-Q12)
# x13 = log(1-Q13)
# x23 = log(1-Q23)
# x123 = log(1-Q123)

# y1 = log(1-P1)
# y2 = log(1-P2)
# y3 = log(1-P3)
# y12 = log(1-P12)
# y13 = log(1-P13)
# y23 = log(1-P23)
# y123 = log(1-P123)

# [y1,y2,y3,y12,y13,y23,y123] = F(x1,x2,x3,x12,x13,x23,x123)

# # Assess the correlation between capacity
# cor(depCap_1, depCap_2) 
# cor(depCap_1, depCap_3) 
# cor(depCap_2, depCap_3) 
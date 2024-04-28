## use this script to generate CPTs for independent cases
library('truncnorm')

#------------------------ Initialize the Working Path ------------------------#
mydirectory <- 'C:\\Users\\taota\\Desktop\\PSAM_seismic\\BayesianNetwork\\BNLearn-independent'

setwd(mydirectory)   # change to mydirectory
options(digits = 12)

#---------------------------------------------------------------------------
## Define a function to generate list of bins
gInt <- function(vX){
  vX_interval = replicate((length(vX)), list()) 
  for (i in 1:(length(vX)-1)) {
    vX_interval[[i]] <- list(smin = vX[i], smax=vX[i+1] )
  }
  return(vX_interval) # component states failure or survive
}

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

#---------------------------------------------------------------------------
#Define component median capacity
median_cap_A = 0.9
median_cap_B = 1.0
median_cap_C = 1.1

#---------------------------------------------------------------------------
#discretize PGA curve
v_gm = c(0, 0.05, 0.25, 0.45, 0.65, 0.85, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50) #TZ
num.gm= length(v_gm)-1 # number of intervals
v_gm_interval = gInt(v_gm)
# v_gm_interval[[1]]$smin
# v_gm_interval[[1]]$smax

#---------------------------------------------------------------------------
##discretize according to normal distribution 
v_sa = c(-2.55, -2.09, -1.62, -1.16, -0.70, -0.23, 0.23, 0.70, 1.16, 1.62, 2.09, 2.55)
v_sb = c(-2.13, -1.74, -1.36, -0.97, -0.58, -0.19, 0.19, 0.58, 0.97, 1.36, 1.74, 2.13)
v_sc = c(-2.16, -1.77, -1.37, -0.98, -0.59, -0.20, 0.20, 0.59, 0.98, 1.37, 1.77, 2.16)


num.sa= length(v_sa)-1 
v_sa_interval = gInt(v_sa)

num.sb= length(v_sb)-1 
v_sb_interval = gInt(v_sb)

num.sc= length(v_sc)-1 
v_sc_interval = gInt(v_sc)

#---------------------------------------------------------------------------
# run MC simulation
#specify number of samples for each replication
samp_size = 1e3

# Specify the CPT for component A
cpt_A = c() 

for (i_sa in 1:num.sa) {
  sa_lower = v_sa_interval[[i_sa]]$smin
  sa_upper = v_sa_interval[[i_sa]]$smax
  samp_sa = rtnorm(samp_size, mean = 0, sd = 0.85,lower=sa_lower, upper=sa_upper)
  
  for (i_GM in 1:num.gm) {
    #sample ground motion
    tmin = v_gm_interval[[i_GM]]$smin
    tmax = v_gm_interval[[i_GM]]$smax

    ########################################################################
    ## Execute Importance Sampling
    mean.proposal = tmax
    sd.proposal = (tmax-tmin)/1
    # Generate ground motion intensity from proposal distribution
    samp_gm <- rtnorm(samp_size,mean.proposal,sd.proposal,lower=tmin, upper=tmax)
    # Determine the weight vector for samples
    W <- p(samp_gm)/q(samp_gm,mean.proposal,sd.proposal,tmin,tmax)
    
    if (i_GM==1) {
      cpt_A = c(cpt_A, 0, 1)
    } else {
      #calculate limit state function
      limit_state_A = exp(log(median_cap_A) + samp_sa) - samp_gm
      # comp_fail_A = sum(limit_state_A<=0)/samp_size #failure
      comp_fail_A = sum(W*as.numeric(limit_state_A<=0))/sum(W)
      comp_suc_A =  1 - comp_fail_A
      cpt_A = c(cpt_A, comp_fail_A, comp_suc_A)
    }
    
  }
  
}

# Specify the CPT for component B
cpt_B = c() 

for (i_sb in 1:num.sb) {
  sb_lower = v_sb_interval[[i_sb]]$smin
  sb_upper = v_sb_interval[[i_sb]]$smax
  samp_sb = rtnorm(samp_size, mean = 0, sd = 0.71,lower=sb_lower, upper=sb_upper)
  
  for (i_GM in 1:num.gm) {
    #sample ground motion
    tmin = v_gm_interval[[i_GM]]$smin
    tmax = v_gm_interval[[i_GM]]$smax

    ########################################################################
    ## Execute Importance Sampling
    #mean.proposal = sqrt(tmin*tmax)
    mean.proposal = tmax
    sd.proposal = (tmax-tmin)/1
    # Generate ground motion intensity from proposal distribution
    samp_gm <- rtnorm(samp_size,mean.proposal,sd.proposal,lower=tmin, upper=tmax)
    # Determine the weight vector for samples
    W <- p(samp_gm)/q(samp_gm,mean.proposal,sd.proposal,tmin,tmax)

    if (i_GM==1) {
      cpt_B = c(cpt_B, 0, 1)
    } else {
      #calculate limit state function
      limit_state_B = exp(log(median_cap_B) + samp_sb) - samp_gm
      # comp_fail_B = sum(limit_state_B<=0)/samp_size #failure
      comp_fail_B = sum(W*as.numeric(limit_state_B<=0))/sum(W)
      comp_suc_B =  1 - comp_fail_B
      
      cpt_B = c(cpt_B, comp_fail_B, comp_suc_B)
    }
    
  }
  
}

# Specify the CPT for component C
cpt_C = c() 

for (i_sc in 1:num.sc) {
  sc_lower = v_sc_interval[[i_sc]]$smin
  sc_upper = v_sc_interval[[i_sc]]$smax
  samp_sc = rtnorm(samp_size, mean = 0, sd = 0.72,lower=sc_lower, upper=sc_upper)
  
  for (i_GM in 1:num.gm) {
    #sample ground motion
    tmin = v_gm_interval[[i_GM]]$smin
    tmax = v_gm_interval[[i_GM]]$smax
    
    ########################################################################
    ## Execute Importance Sampling
    mean.proposal = tmax
    sd.proposal = (tmax-tmin)/1
    # Generate ground motion intensity from proposal distribution
    samp_gm <- rtnorm(samp_size,mean.proposal,sd.proposal,lower=tmin, upper=tmax)
    # Determine the weight vector for samples
    W <- p(samp_gm)/q(samp_gm,mean.proposal,sd.proposal,tmin,tmax)

    if (i_GM==1) {
      cpt_C = c(cpt_C, 0, 1)
    } else {
      #calculate limit state function
      limit_state_C = exp(log(median_cap_C) + samp_sc) - samp_gm
      # comp_fail_C = sum(limit_state_C<=0)/samp_size #failure
      comp_fail_C = sum(W*as.numeric(limit_state_C<=0))/sum(W)
      comp_suc_C =  1 - comp_fail_C
      
      cpt_C = c(cpt_C, comp_fail_C, comp_suc_C)
    }
    
  }
  
}

## output CPTs for independent cases
XX = matrix(cpt_A,nrow=2) 
write.table(XX,"cpt_A_bnlearn_ind.csv", row.names = FALSE, col.names = FALSE, sep=',')

XX = matrix(cpt_B,nrow=2) 
write.table(XX,"cpt_B_bnlearn_ind.csv", row.names = FALSE, col.names = FALSE, sep=',')

XX = matrix(cpt_C,nrow=2) 
write.table(XX,"cpt_C_bnlearn_ind.csv", row.names = FALSE, col.names = FALSE, sep=',')
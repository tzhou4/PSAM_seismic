## use this script to generate CPTs for dependent cases
library('truncnorm')

#------------------------ Initialize the Working Path ------------------------#
mydirectory <- 'C:\\Users\\taota\\Desktop\\PSAM_seismic\\BayesianNetwork\\BNLearn-Dependent'

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
v_sab = c(-1.26, -1.03, -0.80, -0.57, -0.34, -0.11, 0.11, 0.34, 0.57, 0.80, 1.03, 1.26)
v_sac = c(-1.59, -1.30, -1.01, -0.72, -0.43, -0.14, 0.14, 0.43, 0.72, 1.01, 1.30, 1.59)
# v_sbc = 0

num.sab= length(v_sab)-1 
v_sab_interval = gInt(v_sab)

num.sac= length(v_sac)-1 
v_sac_interval = gInt(v_sac)
#---------------------------------------------------------------------------
##discretize according to normal distribution 
v_sa = c(-1.14, -0.93, -0.73, -0.52, -0.31, -0.10, 0.10, 0.31, 0.52, 0.73, 0.93, 1.14)
v_sb = c(-1.50, -1.23, -0.95, -0.68, -0.41, -0.14, 0.14, 0.41, 0.68, 0.95, 1.23, 1.50)
v_sc = c(-1.32, -1.08, -0.84, -0.60, -0.36, -0.12, 0.12, 0.36, 0.60, 0.84, 1.08, 1.32)

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
  samp_sa = rtnorm(samp_size, mean = 0, sd = 0.34,lower=sa_lower, upper=sa_upper)
  
  for (i_sac in 1:num.sac) {
    sac_lower = v_sac_interval[[i_sac]]$smin
    sac_upper = v_sac_interval[[i_sac]]$smax
    samp_sac = rtnorm(samp_size, mean = 0, sd = 0.16,lower=sac_lower, upper=sac_upper)
    
    for (i_sab in 1:num.sab) {
      sab_lower = v_sab_interval[[i_sab]]$smin
      sab_upper = v_sab_interval[[i_sab]]$smax
      samp_sab = rtnorm(samp_size, mean = 0, sd = 0.14,lower=sab_lower, upper=sab_upper)
      
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
          cpt_A = c(cpt_A, 0, 1)
        } else {
          #calculate limit state function
          limit_state_A = exp(log(median_cap_A) + samp_sa + samp_sab + samp_sac) - samp_gm
          #comp_fail_A = sum(limit_state_A<=0)/samp_size #failure
          comp_fail_A = sum(W*as.numeric(limit_state_A<=0))/sum(W)
          comp_suc_A =  1 - comp_fail_A
          cpt_A = c(cpt_A, comp_fail_A, comp_suc_A)
        }
        
      }
      
    }
    
  }
  
}

# Specify the CPT for component B
cpt_B = c() 

for (i_sb in 1:num.sb) {
  sb_lower = v_sb_interval[[i_sb]]$smin
  sb_upper = v_sb_interval[[i_sb]]$smax
  samp_sb = rtnorm(samp_size, mean = 0, sd = 0.29,lower=sb_lower, upper=sb_upper)
  
  for (i_sab in 1:num.sab) {
    sab_lower = v_sab_interval[[i_sab]]$smin
    sab_upper = v_sab_interval[[i_sab]]$smax
    samp_sab = rtnorm(samp_size, mean = 0, sd = 0.14,lower=sab_lower, upper=sab_upper)
    
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
        cpt_B = c(cpt_B, 0, 1)
      } else {
        #calculate limit state function
        limit_state_B = exp(log(median_cap_B) + samp_sb + samp_sab) - samp_gm
        # comp_fail_B = sum(limit_state_B<=0)/samp_size #failure
        comp_fail_B = sum(W*as.numeric(limit_state_B<=0))/sum(W)
        comp_suc_B =  1 - comp_fail_B
        
        cpt_B = c(cpt_B, comp_fail_B, comp_suc_B)
      }
      
    }
    
  }
  
}

# Specify the CPT for component C
cpt_C = c() 

for (i_sc in 1:num.sc) {
  sc_lower = v_sc_interval[[i_sc]]$smin
  sc_upper = v_sc_interval[[i_sc]]$smax
  samp_sc = rtnorm(samp_size, mean = 0, sd = 0.35,lower=sc_lower, upper=sc_upper)
  
  for (i_sac in 1:num.sac) {
    sac_lower = v_sac_interval[[i_sac]]$smin
    sac_upper = v_sac_interval[[i_sac]]$smax
    samp_sac = rtnorm(samp_size, mean = 0, sd = 0.16,lower=sac_lower, upper=sac_upper)
    
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
        limit_state_C = exp(log(median_cap_C) + samp_sc + samp_sac) - samp_gm
        # comp_fail_C = sum(limit_state_C<=0)/samp_size #failure
        comp_fail_C = sum(W*as.numeric(limit_state_C<=0))/sum(W)
        comp_suc_C =  1 - comp_fail_C
        
        cpt_C = c(cpt_C, comp_fail_C, comp_suc_C)
      }
      
    }
    
  }
  
}

## output CPTs
XX = matrix(cpt_A,nrow=2) 
write.table(XX,"cpt_A_bnlearn.csv", row.names = FALSE, col.names = FALSE, sep=',')

XX = matrix(cpt_B,nrow=2) 
write.table(XX,"cpt_B_bnlearn.csv", row.names = FALSE, col.names = FALSE, sep=',')

XX = matrix(cpt_C,nrow=2) 
write.table(XX,"cpt_C_bnlearn.csv", row.names = FALSE, col.names = FALSE, sep=',')

############################################################
## Implement the Bayesian Network Method (Parallel dependent case)
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
library(bnlearn)

#------------------------ Initialize the Working Path ------------------------#
mydirectory <- 'C:\\Users\\taota\\Desktop\\PSAM_seismic\\BayesianNetwork\\BNLearn-Dependent'

setwd(mydirectory)   # change to mydirectory
options(digits = 12)

#--------------------------------------------------------------
#define node labels
#--------------------------------------------------------------
LETTERS = c("EES", "GM", "sAB", "sAC", "sA", "sB", "sC", "A", "B", "C", "S")
#Creating an empty network
seismic = empty.graph(nodes = LETTERS)
class(seismic)

#--------------------------------------------------------------
#Creating a network structure
#--------------------------------------------------------------
# With a specific arc set
arc.set = matrix(c("EES", "GM", 
                   "GM", "A",
                   "GM", "B",
                   "GM", "C",
                   "sAB", "A",
                   "sAB", "B",
                   "sAC", "A",
                   "sAC", "C",
                   "sA", "A",
                   "sB", "B",
                   "sC", "C",
                   "A","S",
                   "B","S",
                   "C","S"),
                 ncol = 2, byrow = TRUE, dimnames = list(NULL, c("from", "to")))
arc.set
#set arcs
arcs(seismic) = arc.set
#plot graph
graphviz.plot(seismic, shape = "ellipse")

#--------------------------------------------------------------
#define conditional probability table
#--------------------------------------------------------------
EES.lv = c("no", "yes")
GM.lv = c("pga0", "pga1", "pga2","pga3","pga4","pga5","pga6","pga7","pga8","pga9","pga10")
sAB.lv = c("sab1", "sab2", "sab3", "sab4", "sab5", "sab6", "sab7", "sab8", "sab9", "sab10", "sab11")
sAC.lv = c("sac1", "sac2", "sac3", "sac4", "sac5", "sac6", "sac7", "sac8", "sac9", "sac10", "sac11")
sA.lv = c("sa1", "sa2", "sa3", "sa4", "sa5", "sa6", "sa7", "sa8", "sa9", "sa10", "sa11")
sB.lv = c("sb1", "sb2", "sb3", "sb4", "sb5", "sb6", "sb7", "sb8", "sb9", "sb10", "sb11")
sC.lv = c("sc1", "sc2", "sc3", "sc4", "sc5", "sc6", "sc7", "sc8", "sc9", "sc10", "sc11")
A.lv = c("Failure", "Success")
B.lv = c("Failure", "Success")
C.lv = c("Failure", "Success")
S.lv = c("Failure", "Success")

#--------------------------------------------------------------
##discretize PGA hazard curve
p_GM_EESno = c(1,0,0,0,0,0,0,0,0,0,0)
p_GM_EESyes = c(0, 9.25E-01, 4.58E-02, 1.30E-02, 5.64E-03, 2.41E-03, 1.14E-03, 9.02E-04, 7.25E-04, 5.93E-04, 4.68E-03)

##discretize according to normal distribution 
p_sab = c(7.05e-3, 2.11e-2, 5.82e-2, 1.20e-1, 1.86e-1, 2.15e-1, 1.86e-1, 1.20e-1, 5.82e-2, 2.11e-2, 7.05e-3)
p_sac = c(7.05e-3, 2.11e-2, 5.82e-2, 1.20e-1, 1.86e-1, 2.15e-1, 1.86e-1, 1.20e-1, 5.82e-2, 2.11e-2, 7.05e-3)

##discretize according to normal distribution 
p_sa = c(7.05e-3, 2.11e-2, 5.82e-2, 1.20e-1, 1.86e-1, 2.15e-1, 1.86e-1, 1.20e-1, 5.82e-2, 2.11e-2, 7.05e-3)
p_sb = c(7.05e-3, 2.11e-2, 5.82e-2, 1.20e-1, 1.86e-1, 2.15e-1, 1.86e-1, 1.20e-1, 5.82e-2, 2.11e-2, 7.05e-3)
p_sc = c(7.05e-3, 2.11e-2, 5.82e-2, 1.20e-1, 1.86e-1, 2.15e-1, 1.86e-1, 1.20e-1, 5.82e-2, 2.11e-2, 7.05e-3)

##run MC simulation to derive the CPT for components A, B and C.
cpt_A <- read.csv(file = 'cpt_A_bnlearn.csv', header = FALSE)
cpt_B <- read.csv(file = 'cpt_B_bnlearn.csv', header = FALSE)
cpt_C <- read.csv(file = 'cpt_C_bnlearn.csv', header = FALSE)

p_a = as.vector(unlist(cpt_A))
p_b = as.vector(unlist(cpt_B))
p_c = as.vector(unlist(cpt_C))

##define CPT for system S
#parallel system
p_s = c(1,0,0,1,0,1,0,1,0,1,0,1,0,1,0,1)

EES.prob = array(c(0.9988, 0.0012), dim = 2, dimnames = list(EES = EES.lv))
GM.prob = array(c(p_GM_EESno,p_GM_EESyes), dim = c(11,2), dimnames = list(GM = GM.lv, EES = EES.lv))

sAB.prob = array(c(p_sab), dim = 11, dimnames = list(sAB = sAB.lv))
sAC.prob = array(c(p_sac), dim = 11, dimnames = list(sAC = sAC.lv))

sA.prob = array(c(p_sa), dim = 11, dimnames = list(sA = sA.lv))
sB.prob = array(c(p_sb), dim = 11, dimnames = list(sB = sB.lv))
sC.prob = array(c(p_sc), dim = 11, dimnames = list(sC = sC.lv))

A.prob = array(c(p_a), 
               dim = c(2,11,11,11,11), dimnames = list(A = A.lv, GM = GM.lv, sAB = sAB.lv,sAC = sAC.lv, sA = sA.lv))
B.prob = array(c(p_b), 
               dim = c(2,11,11,11), dimnames = list(B = B.lv, GM = GM.lv, sAB = sAB.lv, sB = sB.lv))
C.prob = array(c(p_c), 
               dim = c(2,11,11,11), dimnames = list(C = C.lv, GM = GM.lv, sAC = sAC.lv, sC = sC.lv))

S.prob = array(c(p_s), 
               dim = c(2,2,2,2), dimnames = list(S = S.lv, A = A.lv, B = B.lv, C = C.lv))

#SET conditional probability table
cpt = list(EES = EES.prob, GM = GM.prob, 
           sAB= sAB.prob, sAC = sAC.prob, 
           sA = sA.prob, sB = sB.prob, sC = sC.prob,
           A = A.prob, B = B.prob, C = C.prob,
           S = S.prob)

seismic = custom.fit(seismic, cpt)

#--------------------------------------------------------------
# Finally, starting from a bn.fit object we can:
#--------------------------------------------------------------

# generate random samples;
rbn(seismic, n = 10)

#------------------------ Load library ------------------------#
library(bnlearn)
library(gRain)

#------------------------ Inferences ------------------------#
gr.seismic <- as.grain(seismic)
gr.seismic$cptlist$EES
gr.seismic$cptlist$GM

# Compile first
gr.alarm <- compile(gr.seismic) 
gr.alarm

# Query with evidence
gmList = c("pga0", "pga1", "pga2","pga3","pga4","pga5","pga6","pga7","pga8","pga9","pga10")

seismicRisk = c() 
for (gmIndex in 2:11){
  ## prepare evidence
  gmState = gmList[gmIndex]
  ## set evidence
  gr.alarm <- setEvidence(object=gr.alarm, nodes=c('GM'), states=gmState)
  #gr.alarm
  sprob <- querygrain(gr.alarm, type = "joint", nodes = 'S' )
  ## retract evidence
  gr.alarm <- retractEvidence(gr.alarm, nodes=c('GM'))
  #gr.alarm
  seismicRisk = c(seismicRisk, sprob)
}


finalOUTPUT = matrix(seismicRisk,nrow=2)

########## Plots
## graphics.off()
par(mfrow=c(1,2))
## Plot for intersection
plot(finalOUTPUT[1,],type = "l", col = 'red',lty = 1:5, lwd = 1,xlim = NULL, ylim = NULL,
     xlab = 'Acceleration(g)', ylab = 'Conditional Failure Probability')
points(finalOUTPUT[1,],lwd = 3,col = 'red',pch=21)
title(main = 'Failure Probability(Parallel)', sub = NULL,
      cex.main = 0.8)

## Plot for Union
plot(finalOUTPUT[2,],type = "l", col = 'red',lty = 1:5, lwd = 1,xlim = NULL, ylim = NULL,
     xlab = 'Acceleration(g).', ylab = 'Conditional Failure Probability')
points(finalOUTPUT[2,],lwd = 3,col = 'red',pch=21)
title(main = 'Success Probability(Parallel)', sub = NULL,
      cex.main = 0.8)

## output results
write.table(t(finalOUTPUT),"bnlearn_parallel_dep.csv", row.names = FALSE, col.names = FALSE, sep=',')

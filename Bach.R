## ----setup, include=FALSE---------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
options(knitr.purl.inline = TRUE)


## ----load required packages, message=FALSE----------------------------------------
library(markovchain)
library(readr)
library(igraph)
require(reshape2)
require(ggplot2)


## ----dataset, message=F, warning=F------------------------------------------------
jsbach_chorals_harmony <- read_csv("jsbach_chorals_harmony.data", 
                                   col_names = FALSE)
BWV_17.7 = jsbach_chorals_harmony[jsbach_chorals_harmony$X1 == '000306b_',]
bwv177 = factor(BWV_17.7$X17)
head(bwv177)


## ----TM1--------------------------------------------------------------------------
TM = createSequenceMatrix(bwv177, toRowProbs = TRUE); TM


## ----MC1--------------------------------------------------------------------------
MC = as(TM, "markovchain")
#verify that it is a markovchain object
is(MC, 'markovchain')

plot(MC, edge.arrow.size=0.45, main='BWV 17.7 Markov Chain')


## ----2transitions-----------------------------------------------------------------
#initial state
initial_state = c(0,1,0,0,0,0,0,0,0,0,0)
#2 transitions after
after_2_transitions = initial_state * (MC^2); after_2_transitions


## ----properties1------------------------------------------------------------------
#states
states(MC)
#dimension
dim(MC)


## ----transition probabilities-----------------------------------------------------
#beautiful transition, probably it will occur
transitionProbability(MC, 'A_M', 'E_M')
#ugly transition, probably it will not occur 
transitionProbability(MC, 'A#d', 'E_M')


## ----criar um coral---------------------------------------------------------------
#create a coral with 10 harmonic transitions
markovchainSequence(10, MC, t0="A_M")


## ----steadystates1, message=FALSE-------------------------------------------------
#probability with DTMC: stationary distribution
## when the TM is irreducibile
df_MC = melt(data.frame(steadyStates(MC)))
ggplot(df_MC, aes(variable,value), ylim(c(0,0.3))) + geom_point(color='blue') + geom_text(aes(label=round(value, 4)), size=3, vjust=-2)



## ----classifying states-----------------------------------------------------------
#classifying states
transientStates(MC)
absorbingStates(MC)
#identifying recurrent and transient classes
recurrentClasses(MC)
communicatingClasses(MC)
is.accessible(MC, from = "E_M7",to="F#m")
summary(MC)
period(MC)

#converting to igraph
MC.igraph = as(MC,"igraph")
#finding and formatting the clusters
SCC = clusters(MC.igraph, mode="strong") 
V(MC.igraph)$color = rainbow(SCC$no)[SCC$membership]
#plotting
plot(MC.igraph, mark.groups = split(1:vcount(MC.igraph), SCC$membership), main="Communicating classes - strongly connected components")


## ----first passage----------------------------------------------------------------
#first passage time
firstPassage(MC,state = "A_M",5)


## ----2nd chorale------------------------------------------------------------------

bwv2 = jsbach_chorals_harmony[jsbach_chorals_harmony$X1 == '000408b_',]
bwv2 = bwv2$X17
bwv2 = factor(append(bwv2, c('E_m','C_M', 'A_m', 'D_M', 'G_m', 'Eb_M', 'C_m', 'Eb_M', 'G_m', 'G#_M', 'G#_M')))

TM2 = createSequenceMatrix(bwv2, toRowProbs = TRUE)
MC2 = as(TM2, 'markovchain')
summary(MC2)


## ----graph, echo=F----------------------------------------------------------------
#converting to igraph
MC2.igraph = as(MC2,"igraph")
#finding and formatting the clusters
SCC = clusters(MC2.igraph, mode="strong") 
V(MC2.igraph)$color = rainbow(SCC$no)[SCC$membership]
#plotting
plot(MC2.igraph, mark.groups = split(1:vcount(MC2.igraph), SCC$membership), main="Communicating classes - strongly connected components")


## ----Q R I N, message=F-----------------------------------------------------------
#The following function returns the Q, R, and I matrices
extractMatrices = function(mcObj) {
require(matlab)
mcObj = canonicForm(object = mcObj)
  #get the indices of transient and absorbing
transIdx = which(states(mcObj) %in% transientStates(mcObj)) #transient indexes
absIdx = which(states(mcObj) %in% absorbingStates(mcObj)) #absorbing indexes
  #get the Q, R and I matrices
Q = as.matrix(mcObj@transitionMatrix[transIdx,transIdx]) 
R = as.matrix(mcObj@transitionMatrix[transIdx,absIdx])
I = as.matrix(mcObj@transitionMatrix[absIdx, absIdx])
#get the fundamental matrix
N = solve(eye(size(Q)) - Q)
#computing final absortion probabilities 
NR = N %*% R
#return
out = list(
    canonicalForm = mcObj,
    Q = Q,
    R = R,
    I = I,
    N = N,
    NR = NR)
  return(out)
}

#decompose the matrix
MC2_dec = extractMatrices(mcObj = MC2)

#showing the fundamental matrix
MC2_fund = MC2_dec$N; MC2_fund


## ---------------------------------------------------------------------------------
#expected number of steps before being absorbed 
MC2_fund%*%rep(1,17)


## ---------------------------------------------------------------------------------
#calculating B matrix
#the probability to being absorbed in G#_M state as a function of the starting transient state
MC2_B = MC2_fund%*%MC2_dec$R; MC2_B


## ---------------------------------------------------------------------------------
#calculating H, probability of visiting transient state j starting in transient state i
MC2_H = (MC2_fund - matlab::eye(ncol(MC2_fund))) * solve(diag(diag(MC2_fund))); MC2_H


## ---------------------------------------------------------------------------------
period(MC2)


## ---------------------------------------------------------------------------------


steadyStates(MC2)


## ---- message=F-------------------------------------------------------------------
#create R matrix
R = ifelse(TM2<1000,1/dim(MC2))
#apply damping factor
dampened = 0.85*TM2+(1-0.85)*R
dampened = as(dampened, 'markovchain')
df_damp = data.frame(steadyStates(dampened))
df_damp = melt(df_damp)
df_damp$SS = c('dampened')

#stationary distribution
steadyStates(dampened)

#without modification
bwv2_wo = jsbach_chorals_harmony[jsbach_chorals_harmony$X1 == '000408b_',]
bwv2_wo = bwv2_wo$X17
TM2_wo = createSequenceMatrix(bwv2_wo, toRowProbs = TRUE)
MC2_wo = as(TM2_wo, 'markovchain')
original_df = melt(data.frame(steadyStates(MC2_wo)))
original_df$SS = c('original')

#bind the two
df_all = rbind(original_df, df_damp)

#plot
ggplot(df_all, aes(variable,value, color=SS), ylim(c(0,0.15))) + geom_point()



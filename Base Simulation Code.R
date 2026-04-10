##### DATA SIMULATION STUDY - RESEARCH INTERNSHIP #####
# Project: From Symptoms to Systems: CBT Within the Symptom Network
# Student: Sebastian Krebs
# Supervisor: Marie Deserno




################################################################################
############################### Simulation Study ###############################
################################################################################


# Goal: Simulate Expected Disordered Network Based on Theory & Hypotheses
# Goal: Simulate Change In Disordered Network Based on Treatment


################################################################################


################################################################################
#################### Psychopathological Network Simulation #####################
# Type: Theory Involved



############                   Define Functions                    #############

#### parameter (thresholds and weights) function
ind_par <- function(symptoms,
                    auto = FALSE, # define whether symptoms autocorrelate
                    weight_data, # define data to inspire weight generation
                    thres_data,   # define data to inspire threshold generation
                    weight_jitter = 0.05, # define how much weight should 
                    # deviate from data
                    thres_jitter = 0.05,  # define how much thresholds should 
                    connectivity = 0){ 
  # deviate from data
  
  
  
  ### defaults
  ## weight default extracted from Cramer, 2016
  if(missing(weight_data)){
    weight_data <- c(2.1407, 0.7232, 0.1766, 0.2041, 0.2811, 0.6082, 1.129, 0.5763, 2.9840, -0.5389, 0.5217, 0.2392, 0.0000, 3.1650, 0.0000, 1.0530, 0.4273, 0.2045, 0.2672, 0.7033, 0.0000, 0.9409, 0.5311, 0.0000, 0.2041, 0.4020, 0.5475, -0.5009, 0.7484, 0.4459, 0.0000, 0.0000, 0.4724, 0.4890, 1.2951, 0.0000, 0.6849, 0.6564, 0.0000, 0.0000, 0.2219, 0.2914, 0.0000, 0.4048, 0.0000, 1.0979, 0.5070, 0.0000, 0.4112, 0.2284, 0.4546, 0.8279, 1.4768, 0.3751, 1.5718, 1.8733, 0.6826, 0.3772, 0.5226, 0.1203, 0.0000, 0.0000, 0.2708, 0.3893, 0.3491, 0.2362, 1.0211, 0.8178, 0.1063, 0.0000, 0.4177, 0.0000, 0.2585, 0.0597, 0.9414, 0.7233, 0.4935, 0.6660, 2.0693, 0.4986, 0.2520, 0.0000, 0.1198, 0.1620, 0.4514, 0.2151, 0.1939, 0.1407, 0.0000, 1.4769, 0.2156)
  }
  
  ## threshold default extracted from Cramer, 2016
  if(missing(thres_data)){
    thres_data <- c(-3.1946, -4.3092, -3.8332, -3.9153, -3.9012, -3.0246, -4.4480, -3.1753, -4.3372, -2.8269, -4.4272, -4.0421, -5.8303)
  }
  
  
  
  ### warnings and checks
  if(!is.vector(symptoms) || !is.character(symptoms)){
    stop("symptoms should be reported in a single character vector")
  }
  if(auto == TRUE){
    warning("you are currently assuming that symptoms influence themselves over time")
  }
  if(!is.vector(weight_data) || !is.numeric(weight_data)){
    stop("weigth_data should be reported in a single numeric vector")
  }
  if(!is.vector(thres_data) || !is.numeric(thres_data)){
    stop("weigth_data should be reported in a single numeric vector")
  }
  if(!is.numeric(weight_jitter) || !is.numeric(thres_jitter)){
    stop("weight_jitter and thres_jitter should be numeric")
  }
  
  
  
  ### define generating functions (see above)
  ## define function to simulate new weights based on existing data
  weight_sim <- function(data, 
                         n = length(data), 
                         jitter_sd = 0.05, 
                         connectivity) {
    
    # checks
    if(!is.vector(data)){
      stop("data needs to be a prepared vector excluding the 0 diagonal and including corresponding weights only once")
    }
    
    # resample from data
    resampled <- sample(data, size = n, replace = TRUE)
    
    
    # select only nonzeros (zeros should be maintained, as absence of edge)
    nonzero <- which(resampled != 0)
    
    # add connectivity
    increase <- rep(0, length(nonzero))
    
    if(connectivity != 0){
      increase <- rnorm(length(nonzero), connectivity)
    }
    
    # add jitter only if there are non-zero values 
    if (length(nonzero) > 0) {
      noise <- rnorm(length(nonzero), mean = 0, sd = jitter_sd)
      resampled[nonzero] <- resampled[nonzero] + noise + increase
      
      
      # jittered data should remain within original data bounds
      resampled[nonzero] <- pmax(resampled[nonzero], min(data))
      resampled[nonzero] <- pmin(resampled[nonzero], max(data))
    }
    
    
    # return new data
    return(resampled)
  }
  
  
  ## define function to simulate thresholds, similar to the weight function
  thres_sim <- function(data, 
                        n = length(data), 
                        jitter_sd = 0.05) {
    
    # checks
    if(!is.vector(data)){
      stop("data needs to be a vector including all datapoints")
    }
    
    # resample from data
    resampled <- sample(data, size = n, replace = TRUE)
    
    # add noise to data
    noise <- rnorm(n, mean = 0, sd = jitter_sd)
    resampled <- resampled + noise
    
    # jittered data should remain within original data bounds
    resampled <- pmax(resampled, min(data))
    resampled <- pmin(resampled, max(data))
    
    # return new data
    return(resampled)
  }
  
  
  
  ### generate unique participant data
  ## generate weight matrix
  # calculate weights
  weight <- weight_sim(data = weight_data, 
                       n = choose(length(symptoms), 2) + length(symptoms), 
                       # generate enough weights to fill matrix above the 
                       # diagonal + the diagonal
                       jitter_sd = weight_jitter,
                       connectivity = connectivity)
  
  # create mirrored matrix filled with weights
  w <- matrix(0, length(symptoms), length(symptoms))
  w[upper.tri(w)] <- weight[1:choose(length(symptoms), 2)]
  weights <- w + t(w)
  
  # if requested set diagonal (autocorrelations)
  if(auto == TRUE){
    diag(weights) <- weight[(choose(length(symptoms), 2) + 1):
                              (choose(length(symptoms), 2) + length(symptoms))]
  }
  
  # assign sympom names
  rownames(weights) <- symptoms
  colnames(weights) <- symptoms
  
  
  ## generate threshold vector
  thresh <- thres_sim(data = thres_data,
                      n = length(symptoms),
                      jitter_sd = thres_jitter)
  
  
  
  ### create individual personalized dataframe
  pers.par <- data.frame(Thresholds = thresh, weights)
  
  
  
  ### return dataframe as output
  return(pers.par)
}




#### total activation function (Cramer, 2016)
A <- function(weights,    # W: vector of weights between current node and other 
              # nodes
              neighbors){ # X: vector of states of all other nodes
  
  
  ### warning and checks
  warning("weights should be aligned with the corresponding neighbors")
  
  if(!is.vector(weights) || !is.numeric(weights)){
    stop("weights should be a numeric vector")
  }
  if(!is.vector(neighbors) || !is.numeric(neighbors)){
    stop("neighbors should be a numeric vector")
  }
  if(length(weights) != length(neighbors)){
    stop("weights and neighbors should align in length")
  }
  
  
  ### function calculation
  A <- sum(weights * neighbors)
  
  
  ### return A
  return(A)
}




#### probability function (Cramer, 2016)
P <- function(thresh, # b, estimated threshold
              act){   # a, estimated activation
  
  ### checks
  if(!is.numeric(thresh) || length(thresh) > 1){
    stop("thresh should be numeric and of length 1")
  }
  if(!is.numeric(act) || length(act) > 1){
    stop("act should be numeric and of length 1")
  }
  
  
  ### absolute value of b
  thresh <- abs(thresh)
  
  
  ### calculate probability of spin
  p <- 1/(1 + exp(thresh - act))
  
  
  ### return spin probability
  return(p)
}




#### base simulation function (Cramer, 2016 + Glauber Dynamics)
ind_sim <- function(pars,         # indicate the personal simulated parameters
                    steps = 1000, # number of MC steps (all symptoms on average)
                    days,         # alternative to steps
                    spdays,       # steps per days; default is number of nodes
                    init = NULL){
  
  ### checks
  if(!is.data.frame(pars)){
    stop("pars should be a dataframe including only thresholds and weights")
  }
  warning("make sure thresholds are in the first column and weights in the following")
  
  if(!exists("A") || !is.function(A)){
    stop("activation function A() must be loaded")
  }
  if(!exists("P") || !is.function(P)){
    stop("probability function P() must be loaded")
  }
  
  
  ### extract parameters
  weights <- as.matrix(pars[, -1])
  thresh <- pars[, 1]
  
  n_symptoms <- nrow(weights)
  
  
  ## extract steps
  if(missing(days) || missing(spdays)){
    
    # check
    if(!is.numeric(steps) || length(steps) > 1){
      stop("steps should be a numeric object of lenght 1")
    }
    
    # set steps
    steps <- steps
    
    # calculate steps if indicated via days
  } else {
    
    # check
    if(!is.numeric(days) || length(days) > 1){
      stop("days should be a numeric object of lenght 1")
    }
    if(!is.numeric(spdays) || length(spdays) > 1){
      stop("spdays should be a numeric object of lenght 1")
    }
    
    # set steps
    steps <- days*spdays
  }
  
  
  ## initial state
  if(is.null(init)){
    
    # set state
    state <- rbinom(n_symptoms, 1, 0.5)         # is 0.5 sensible here?
    
  } else {
    
    # check
    if(!is.numeric(init) || !is.vector(init)){
      stop("init should be a numeric vector")
    }
    
    # set state
    state <- init
  }
  
  
  ## initiate matrix collecting states
  states <- matrix(NA, steps, n_symptoms)
  
  
  
  ### updating of state using A and P functions and Glauber dynamics (MC steps)
  for(t in 1:steps){
    
    
    ## calculate changes for one MC step (on average all symptoms)
    for (u in 1:n_symptoms){
      
      # choose node randomly
      i <- sample(1:n_symptoms, 1)       # should I include updating probabilities as vicious circles
      
      # activation of node
      act <- A(weights[i, ], state)
      
      # probability of spin
      p <- P(thresh[i], act)
      
      # update node
      state[i] <- rbinom(1, 1, p)
    }
    
    
    ## updated states
    states[t, ] <- state
  }
  
  
  
  ### add appropriate symptom names
  colnames(states) <- rownames(pars)
  
  
  
  ### return states
  return(states)
}





############                  Disordered Network                   #############




## Disordered Network



## Intervened Network LBA


## Intervened Network KBA


## Intervened Network CT





################################################################################
############################## Simulation Analysis #############################
################################################################################


# Goal: Check Whether Expected Dyamics are Found
# Goal: Check What Dynamics Are Found at the Cross-Sectional Level


################################################################################
### Within Person Analysis ###


### Construct Cross-Sectional Networks ###


### Network Intervention Analysis ###


### Network Comparison Test ###

## motif analyses ##

## density measures ##




################################################################################
############## Data Analysis via Comparison to Simulated Networks ##############
################################################################################



### Construct Cross-Sectional Networks ###

### Network Intervention Analysis ###

### Network Comparison Test ###






################################################################################
# NOTES


### Trial Self-Built Ising Model ###

# 'weights' is the network structure 



Hamiltonian <- function(spin,
                        weights,
                        field = 0){
  
  n <- length(spin)
  
  # external field
  f_inf <- sum(field * spin)
  
  # interactions
  w_inf <- 0
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      w_inf <- w_inf + weights[i,j]*spin[i]*spin[j]
    }
  }
  
  H <- -f_inf - w_inf
  return(H)
}

p_conf <- function(temp = 1, H, spins, weights, field){
  
  X <- as.matrix(expand.grid(rep(list(c(0,1)), length(spins))))
  n_conf <- nrow(X)
  
  Z <- 0
  
  for(i in 1:n_conf){
    Z <- Z + exp(-temp * Hamiltonian(X[i,], weights, field))
  }
  
  p <- exp(-temp*H)/Z
  return(p)
}

spin <- c(0, 1, 1, 0, 0, 0, 1, 1, 0)
set.seed(1)

weights <- matrix(runif(9*9, -1, 1), nrow = 9)

# remove self-interactions
diag(weights) <- 0

# make symmetric (typical Ising assumption)
weights <- (weights + t(weights))/2


field <- runif(9, -0.5, 0.5)

H <- Hamiltonian(
  spin = spin,
  weights = weights,
  field = field
)

H

p_conf(
  temp = 1,
  H = H,
  spins = spin,
  weights = weights,
  field = field
)





library(tidyverse)

set.seed(123)

## create a small world network
n_nodes <- 9
sw_network <- matrix(rbinom(n_nodes^2,1,0.3), n_nodes, n_nodes)
sw_network[lower.tri(sw_network)] <- t(sw_network)[lower.tri(sw_network)] # symmetric

## define parameters
alpha <- runif(n_nodes, -0.5, 0.5)  # external fields
beta_seq <- seq(0, 1, 0.05)         # sweeping beta(temp)
n_sample <- 1000
n_iter <- 100                       # Gibbs iterations per sample

## Hamiltonian
Hamiltonian <- function(spin,
                        weights,
                        field = 0){
  
  n <- length(spin)
  
  # external field
  f_inf <- sum(field * spin)
  
  # interactions
  w_inf <- 0
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      w_inf <- w_inf + weights[i,j]*spin[i]*spin[j]
    }
  }
  
  H <- -f_inf - w_inf
  return(H)
}


## Gibbs sampler for one step 
gibbs_step <- function(spin, weights, field, beta){
  n <- length(spin)
  for(i in 1:n){
    
    # compute energy difference if we flip spin i
    spin_flip <- spin
    spin_flip[i] <- 1 - spin_flip[i]  # assume spins are 0/1
    dH <- Hamiltonian(spin_flip, weights, field) - Hamiltonian(spin, weights, field)
    p_flip <- 1 / (1 + exp(beta * dH)) # logistic probability
    if(runif(1) < p_flip){
      spin[i] <- spin_flip[i]
    }
  }
  return(spin)
}


## run simulation
mean_spins <- data.frame()

for(beta in beta_seq){
  
  # initialize spins randomly
  spin <- sample(c(0,1), n_nodes, replace = TRUE)
  spins_mat <- matrix(NA, nrow = n_sample, ncol = n_nodes)
  
  for(s in 1:n_sample){
    for(iter in 1:n_iter){
      spin <- gibbs_step(spin, sw_network, alpha, beta)
    }
    spins_mat[s, ] <- spin
  }
  
  mean_spin <- rowMeans(spins_mat)
  mean_spins <- rbind(mean_spins, data.frame(beta = beta, mean_spin = mean_spin))
}


## plot "bifurcation"
mean_spins %>%
  ggplot(aes(x = beta, y = mean_spin)) +
  geom_jitter(width = 0.01, alpha = 0.3, color = "blue") +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  theme_classic() +
  ggtitle("Pseudo-bifurcation in a 9-node Ising network (custom sampler)") +
  xlab(expression(beta)) +
  ylab("Mean Spin")

### Finneman, Theoretical and Statistical Ising Model
#install.packages("IsingSampler")
#install.packages("bayess")
#install.packages("igraph")
#install.packages("parSim")
#install.packages("ggplot2")

library(IsingSampler)
library(bayess)
library(igraph)
library(parSim)
library(ggplot2)

















library(IsingSampler)
library(tidyverse)

set.seed(123)

# 9-node small-world network example
n_nodes <- 9
sw_network <- matrix(rbinom(n_nodes^2,1,0.3), n_nodes, n_nodes)
sw_network[lower.tri(sw_network)] <- t(sw_network)[lower.tri(sw_network)] # symmetric

# Parameters
alpha <- 0   # external field
beta_seq <- seq(0, 1, 0.05)  # sweeping beta
n_sample <- 1000

# Storage for mean spins
mean_spins <- data.frame()

for(beta in beta_seq){
  samp <- IsingSampler(
    n = n_sample,
    graph = sw_network,
    thresholds = rep(alpha, n_nodes),
    beta = beta,
    nIter = 100,
    responses = c(-1,1)
  )
  
  mean_spin <- apply(samp, 1, mean)
  mean_spins <- rbind(mean_spins,
                      data.frame(beta = beta, mean_spin = mean_spin))
}

# plot “pseudo-bifurcation”
mean_spins %>%
  ggplot(aes(x = beta, y = mean_spin)) +
  geom_jitter(width = 0.01, alpha = 0.3, color = "blue") +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  theme_classic() +
  ggtitle("Pseudo-bifurcation in a 9-node Ising network") +
  xlab(expression(beta)) +
  ylab("Mean Spin")





##### Cramer's (2016) formalisation #####
#### network characteristics
### define weights and thresholds (from empirical data but randomised for indiv.)
## define function to simulate new weights based on existing data
weight_sim <- function(data, n = length(data), jitter_sd = 0.05) {
  
  # checks
  if(!is.vector(data)){
    stop("data needs to be a prepared vector excluding the 0 diagonal and including corresponding weights only once")
  }
  
  # resample from data
  resampled <- sample(data, size = n, replace = TRUE)
  
  
  # select only nonzeros (zeros should be maintained, as absence of edge)
  nonzero <- which(resampled != 0)
  
  
  # add jitter only if there are non-zero values 
  if (length(nonzero) > 0) {
    noise <- rnorm(length(nonzero), mean = 0, sd = jitter_sd)
    resampled[nonzero] <- resampled[nonzero] + noise
    
    
    # jittered data should remain within original data bounds
    resampled[nonzero] <- pmax(resampled[nonzero], min(data))
    resampled[nonzero] <- pmin(resampled[nonzero], max(data))
  }
  return(resampled)
}


## define function to simulate thresholds, similar to the weight function
thres_sim <- function(data, n = length(data), jitter_sd = 0.05) {
  
  # checks
  if(!is.vector(data)){
    stop("data needs to be a vector including all datapoints")
  }
  
  # resample from data
  resampled <- sample(data, size = n, replace = TRUE)
  
  # add noise to data
  noise <- rnorm(n, mean = 0, sd = jitter_sd)
  resampled <- resampled + noise
  
  # jittered data should remain within original data bounds
  resampled <- pmax(resampled, min(data))
  resampled <- pmin(resampled, max(data))
  
  # return new data
  return(resampled)
}



# prepare empirical data to use for weights
w <- read.csv("EmpiricalWeightParametersCramer2016.txt", header = TRUE, sep = "")
w_vec <- w[upper.tri(w)] # store each weight one time (excluding the 0 diagonal)


# prepare empirical data to use for thresholds
thres <- read.csv("EmpiricalThresholdParametersCramer2016.txt", header = TRUE, sep = "")
thres_vec <- thres[, 1]



# simulate new data using function and comparing plots
par(mfrow = c(1, 2))

plot(weight_sim(w_vec))
plot(w_vec)


# simulate new data using function and comparing plots
plot(thres_sim(thres_vec, 10))
plot(thres_vec)



## define symptoms
symptoms <- c("Depressed Mood", 
              "Pleasure/Interest",
              "Weight",
              "Insomnia/Hypersomnia", 
              "Psychomotor Agitation/Retardation",
              "Fatigue",
              "Worthlessness/Guilt",
              "Concentration",
              "Suicidality")



#### define function to create dataframe with characteristics for one person
ind_par <- function(symptoms,
                    auto = FALSE, # define whether symptoms autocorrelate
                    weight_data, # define data to inspire weight generation
                    thres_data,   # define data to inspire threshold generation
                    weight_jitter = 0.05, # define how much weight should 
                    # deviate from data
                    thres_jitter = 0.05){ # define how much thresholds should 
  # deviate from data
  
  
  
  ### defaults
  ## weight default extracted from Cramer, 2016
  if(missing(weight_data)){
    weight_data <- c(2.1407, 0.7232, 0.1766, 0.2041, 0.2811, 0.6082, 1.129, 0.5763, 2.9840, -0.5389, 0.5217, 0.2392, 0.0000, 3.1650, 0.0000, 1.0530, 0.4273, 0.2045, 0.2672, 0.7033, 0.0000, 0.9409, 0.5311, 0.0000, 0.2041, 0.4020, 0.5475, -0.5009, 0.7484, 0.4459, 0.0000, 0.0000, 0.4724, 0.4890, 1.2951, 0.0000, 0.6849, 0.6564, 0.0000, 0.0000, 0.2219, 0.2914, 0.0000, 0.4048, 0.0000, 1.0979, 0.5070, 0.0000, 0.4112, 0.2284, 0.4546, 0.8279, 1.4768, 0.3751, 1.5718, 1.8733, 0.6826, 0.3772, 0.5226, 0.1203, 0.0000, 0.0000, 0.2708, 0.3893, 0.3491, 0.2362, 1.0211, 0.8178, 0.1063, 0.0000, 0.4177, 0.0000, 0.2585, 0.0597, 0.9414, 0.7233, 0.4935, 0.6660, 2.0693, 0.4986, 0.2520, 0.0000, 0.1198, 0.1620, 0.4514, 0.2151, 0.1939, 0.1407, 0.0000, 1.4769, 0.2156)
  }
  
  ## threshold default extracted from Cramer, 2016
  if(missing(thres_data)){
    thres_data <- c(-3.1946, -4.3092, -3.8332, -3.9153, -3.9012, -3.0246, -4.4480, -3.1753, -4.3372, -2.8269, -4.4272, -4.0421, -5.8303)
  }
  
  
  
  ### warnings and checks
  if(!is.vector(symptoms) || !is.character(symptoms)){
    stop("symptoms should be reported in a single character vector")
  }
  if(auto == TRUE){
    warning("you are currently assuming that symptoms influence themselves over time")
  }
  if(!is.vector(weight_data) || !is.numeric(weight_data)){
    stop("weigth_data should be reported in a single numeric vector")
  }
  if(!is.vector(thres_data) || !is.numeric(thres_data)){
    stop("weigth_data should be reported in a single numeric vector")
  }
  if(!is.numeric(weight_jitter) || !is.numeric(thres_jitter)){
    stop("weight_jitter and thres_jitter should be numeric")
  }
  
  
  
  ### define generating functions (see above)
  ## define function to simulate new weights based on existing data
  weight_sim <- function(data, n = length(data), jitter_sd = 0.05) {
    
    # checks
    if(!is.vector(data)){
      stop("data needs to be a prepared vector excluding the 0 diagonal and including corresponding weights only once")
    }
    
    # resample from data
    resampled <- sample(data, size = n, replace = TRUE)
    
    
    # select only nonzeros (zeros should be maintained, as absence of edge)
    nonzero <- which(resampled != 0)
    
    
    # add jitter only if there are non-zero values 
    if (length(nonzero) > 0) {
      noise <- rnorm(length(nonzero), mean = 0, sd = jitter_sd)
      resampled[nonzero] <- resampled[nonzero] + noise
      
      
      # jittered data should remain within original data bounds
      resampled[nonzero] <- pmax(resampled[nonzero], min(data))
      resampled[nonzero] <- pmin(resampled[nonzero], max(data))
    }
    return(resampled)
  }
  
  
  ## define function to simulate thresholds, similar to the weight function
  thres_sim <- function(data, n = length(data), jitter_sd = 0.05) {
    
    # checks
    if(!is.vector(data)){
      stop("data needs to be a vector including all datapoints")
    }
    
    # resample from data
    resampled <- sample(data, size = n, replace = TRUE)
    
    # add noise to data
    noise <- rnorm(n, mean = 0, sd = jitter_sd)
    resampled <- resampled + noise
    
    # jittered data should remain within original data bounds
    resampled <- pmax(resampled, min(data))
    resampled <- pmin(resampled, max(data))
    
    # return new data
    return(resampled)
  }
  
  
  
  ### generate unique participant data
  ## generate weight matrix
  # calculate weights
  weight <- weight_sim(data = weight_data, 
                       n = choose(length(symptoms), 2) + length(symptoms), 
                       # generate enough weights to fill matrix above the 
                       # diagonal + the diagonal
                       jitter_sd = weight_jitter)
  
  # create mirrored matrix filled with weights
  w <- matrix(0, length(symptoms), length(symptoms))
  w[upper.tri(w)] <- weight[1:choose(length(symptoms), 2)]
  weights <- w + t(w)
  
  # if requested set diagonal (autocorrelations)
  if(auto == TRUE){
    diag(weights) <- weight[(choose(length(symptoms), 2) + 1):
                              (choose(length(symptoms), 2) + length(symptoms))]
  }
  
  # assign sympom names
  rownames(weights) <- symptoms
  colnames(weights) <- symptoms
  
  
  ## generate threshold vector
  thresh <- thres_sim(data = thres_data,
                      n = length(symptoms),
                      jitter_sd = thres_jitter)
  
  
  
  ### create individual personalized dataframe
  pers.par <- data.frame(Thresholds = thresh, weights)
  
  
  
  ### return dataframe as output
  return(pers.par)
}



#### total activation function ####


A <- function(weights,    # W: vector of weights between current node and other 
              # nodes
              neighbors){ # X: vector of states of all other nodes
  
  
  ### warning and checks
  warning("weights should be aligned with the corresponding neighbors")
  
  if(!is.vector(weights) || !is.numeric(weights)){
    stop("weights should be a numeric vector")
  }
  if(!is.vector(neighbors) || !is.numeric(neighbors)){
    stop("neighbors should be a numeric vector")
  }
  if(length(weights) != length(neighbors)){
    stop("weights and neighbors should align in length")
  }
  
  
  ### function calculation
  A <- sum(weights * neighbors)
  
  
  ### return A
  return(A)
}



#### probability function ####
P <- function(thresh, # b, estimated threshold
              act){   # a, estimated activation
  
  ### checks
  if(!is.numeric(thresh) || length(thresh) > 1){
    stop("thresh should be numeric and of length 1")
  }
  if(!is.numeric(act) || length(act) > 1){
    stop("act should be numeric and of length 1")
  }
  
  
  ### absolute value of b
  thresh <- abs(thresh)
  
  
  ### calculate probability of spin
  p <- 1/(1 + exp(thresh - act))
  
  
  ### return spin probability
  return(p)
}


pp <- ind_par(symptoms)


#### Simulation of Person ####
ind_sim <- function(pars,         # indicate the personal simulated parameters
                    steps = 1000, # number of MC steps (all symptoms on average)
                    days,         # alternative to steps
                    spdays,       # steps per days; default is number of nodes
                    init = NULL){
  
  ### checks
  if(!is.data.frame(pars)){
    stop("pars should be a dataframe including only thresholds and weights")
  }
  warning("make sure thresholds are in the first column and weights in the following")
  
  if(!exists("A") || !is.function(A)){
    stop("activation function A() must be loaded")
  }
  if(!exists("P") || !is.function(P)){
    stop("probability function P() must be loaded")
  }
  
  
  ### extract parameters
  weights <- as.matrix(pars[, -1])
  thresh <- pars[, 1]
  
  n_symptoms <- nrow(weights)
  
  
  ## extract steps
  if(missing(days) || missing(spdays)){
    
    # check
    if(!is.numeric(steps) || length(steps) > 1){
      stop("steps should be a numeric object of lenght 1")
    }
    
    # set steps
    steps <- steps
    
    # calculate steps if indicated via days
  } else {
    
    # check
    if(!is.numeric(days) || length(days) > 1){
      stop("days should be a numeric object of lenght 1")
    }
    if(!is.numeric(spdays) || length(spdays) > 1){
      stop("spdays should be a numeric object of lenght 1")
    }
    
    # set steps
    steps <- days*spdays
  }
  
  
  ## initial state
  if(is.null(init)){
    
    # set state
    state <- rbinom(n_symptoms, 1, 0.2)         # is 0.2 sensible here?
    
  } else {
    
    # check
    if(!is.numeric(init) || !is.vector(init)){
      stop("init should be a numeric vector")
    }
    
    # set state
    state <- init
  }
  
  
  ## initiate matrix collecting states
  states <- matrix(NA, steps, n_symptoms)
  
  
  
  ### updating of state using A and P functions and Glauber dynamics (MC steps)
  for(t in 1:steps){
    
    
    ## calculate changes for one MC step (on average all symptoms)
    for (u in 1:n_symptoms){
      
      # choose node randomly
      i <- sample(1:n_symptoms, 1)       # should I include updating probabilities as vicious circles
      
      # activation of node
      act <- A(weights[i, ], state)
      
      # probability of spin
      p <- P(thresh[i], act)
      
      # update node
      state[i] <- rbinom(1, 1, p)
    }
    
    
    ## updated states
    states[t, ] <- state
  }
  
  
  
  ### add appropriate symptom names
  colnames(states) <- rownames(pars)
  
  
  
  ### return states
  return(states)
}


p <- ind_sim(pp, days = 60, spdays = 3)
p


# we need asynchronous updating



# how do we implement the asynchronous updating?


# how do we implement feedback loops?






#### PHQ 9 TRANSLATION ####
# Question: Over the last 2 weeks, how often have you been bothered by any of 
#           the following problems?
# Answers:  Not at all (0), Several Days (1), More Than Half the Days (2), 
#           Nearly Every Day (3)




#### add self-loops to represent persistance of symptoms
#### simulate individual time scales per symptom





















#DISORDER NETWORK



#### SetUp


### define nodes
symptoms <- c("Pleasure/Interest",
              "Depressed Mood",
              "Insomnia/Hypersomnia", 
              "Fatigue",
              "Weight",
              "Worthlessness/Guilt",
              "Concentration",
              "Psychomotor Agitation/Retardation",
              "Suicidality")


#### Simulate Data
### initiate participants and dataframe to save data
part <- 900
steps <- 100                  # step is interpreted as one day (opportunity to change for all symptoms)

data <- data.frame()

PHQ9 <- data.frame(matrix(NA,
                          nrow = part,
                          ncol = length(symptoms)))
colnames(PHQ9) <- symptoms


### calculate data per participant over time
for(i in 1:part){
  
  
  ## determine parameter per participant in pars_data
  pars_data <- ind_par(symptoms, auto = TRUE, connectivity = 1)
  
  
  ## simulate data for 1000 steps
  sim_data <- ind_sim(pars_data, steps = steps)
  
  # add names
  rownames(sim_data) <- paste0("Sim", i, " Step", 1:steps)
  
  # add to cumulative dataframe
  data <- rbind(data, sim_data)
}



### Estimate PHQ9 Scores per Participant For the Last 14 Days (Steps)
# select part of df for participants
for(i in 1:part){
  
  part_data <- data[(1 + (i - 1) * steps):(i * steps), ]
  
  
  # select the last 14 steps (days of symptoms) and assign scores
  symp_data <- tail(part_data, 14)
  
  
  for(j in 1:length(symptoms)) {
    
    # assign 0, if symptom is present less than 25% of the time
    if(sum(symp_data[, j]) <= length(symp_data) * 0.25) {
      PHQ9[i, j] <- 0
    } else if (sum(symp_data[, j]) <= length(symp_data) * 0.5) {
      PHQ9[i, j] <- 1
    } else if (sum(symp_data[, j]) <= length(symp_data) * 0.75) {
      PHQ9[i, j] <- 2
    } else {
      PHQ9[i, j] <- 3
    }
  }
  
  # adjust column names
  rownames(PHQ9)[i] <- paste0("Sim", i, " PHQ9")
}

##### DEFINE CB ASSUMPTIONS and CBT PREDICTIONS #####

# Answers:  Not at all (0), Several Days (1), More Than Half the Days (2), 
#           Nearly Every Day (3)

##### DATA SIMULATION STUDY - RESEARCH INTERNSHIP #####
# Project: From Symptoms to Systems: CBT Within the Symptom Network
# Student: Sebastian Krebs
# Supervisor: Marie Deserno


         
         
################################################################################
########################## Maintenance Formalisations ##########################
################################################################################ 

# Explanation:  Each CBT treatment has assumed underlying maintenance mechanisms 
#               of MDD within CBT. Maintenance mechanisms in behavioral 
#               activation by Lewinsohn et al. (1986; LBA), Kanter et al. (2008; 
#               KBA) and cognitive therapy (Beck et al., 1979; CT) are
#               implemented formalised below.

# Formalization:  For now, to maintain simplicity of the models, only thresholds
#               and connectivity will be formalized. Implementing novel 
#               functions does indeed seem appropriate, but is beyond the scope
#               of this project, will improve model complexity significantly and
#               requires large amounts of supporting evidence.




################################################################################
######          formalizable base weight and threshold functions          ######

# Assumption:   symptoms correspond to MDD symptoms as defined under the DSM-5
#               (see Borsboom, 2017)



#### parameter (thresholds and weights) function
ind_par <- function(auto = FALSE,          # symptom autocorrelation
                    
                    weight_data,           # reference weight generation data
                    weight_jitter = 0.05,  # deviation from generation data
                    
                    # include CBT treatment maintenance weights
                    LBA_weights = FALSE,
                    KBA_weights = FALSE,
                    CT_weights = FALSE,
                    
                    thres_data,            # reference threshold generation data
                    thres_jitter = 0.05,   # define how much thresholds should 
                    
                    # include CBT treatment maintenance thresholds
                    LBA_thres = FALSE,
                    KBA_thres = FALSE,
                    CT_thres = FALSE,
                    
                    connectivity = NULL){  # character - low, moderate, or high 
                                           # (Cramer et al., 2016)
  
  
  ### data defaults
  ## weight default (from Cramer et al., 2016)
  if(missing(weight_data)){
    weight_data <- c(2.1407, 0.7232, 0.1766, 0.2041, 0.2811, 0.6082, 1.129, 0.5763, 2.9840, -0.5389, 0.5217, 0.2392, 0.0000, 3.1650, 0.0000, 1.0530, 0.4273, 0.2045, 0.2672, 0.7033, 0.0000, 0.9409, 0.5311, 0.0000, 0.2041, 0.4020, 0.5475, -0.5009, 0.7484, 0.4459, 0.0000, 0.0000, 0.4724, 0.4890, 1.2951, 0.0000, 0.6849, 0.6564, 0.0000, 0.0000, 0.2219, 0.2914, 0.0000, 0.4048, 0.0000, 1.0979, 0.5070, 0.0000, 0.4112, 0.2284, 0.4546, 0.8279, 1.4768, 0.3751, 1.5718, 1.8733, 0.6826, 0.3772, 0.5226, 0.1203, 0.0000, 0.0000, 0.2708, 0.3893, 0.3491, 0.2362, 1.0211, 0.8178, 0.1063, 0.0000, 0.4177, 0.0000, 0.2585, 0.0597, 0.9414, 0.7233, 0.4935, 0.6660, 2.0693, 0.4986, 0.2520, 0.0000, 0.1198, 0.1620, 0.4514, 0.2151, 0.1939, 0.1407, 0.0000, 1.4769, 0.2156)
  }
  
  ## threshold default (from Cramer et al., 2016)
  if(missing(thres_data)){
    thres_data <- c(-3.1946, -4.3092, -3.8332, -3.9153, -3.9012, -3.0246, -4.4480, -3.1753, -4.3372, -2.8269, -4.4272, -4.0421, -5.8303)
  }
  
  
  
  ### warnings and checks
  if(auto == TRUE){
    warning("you are currently assuming that symptoms influence themselves over time")
  }
  
  if(!is.vector(weight_data) || !is.numeric(weight_data)){
    stop("weight_data should be reported in a single numeric vector")
  }
  if(!is.vector(thres_data) || !is.numeric(thres_data)){
    stop("weight_data should be reported in a single numeric vector")
  }
  if(!is.numeric(weight_jitter) || !is.numeric(thres_jitter)){
    stop("weight_jitter and thres_jitter should be numeric")
  }
  
  if(!is.logical(LBA_weights) || !is.logical(KBA_weights) || !is.logical(CT_weights)){
    stop("use logicals for inclusion of maintenance weights")
  }
  if(length(LBA_weights) != 1 || length(KBA_weights) != 1 || length(CT_weights) != 1){
    stop("indicate maintenance weights inclusion with 1(!) logical")
  }
  
  if(!is.logical(LBA_thres) || !is.logical(KBA_thres) || !is.logical(CT_thres)){
    stop("use logicals for inclusion of maintenance thresholds")
  }
  if(length(LBA_thres) != 1 || length(KBA_thres) != 1 || length(CT_thres) != 1){
    stop("indicate maintenance thresholds inclusion with 1(!) logical")
  }

  
  
  ### base functions to simulate parameters
  
  
  ## weight function
  weight_sim <- function(data, 
                         n = length(data), 
                         jitter_sd = 0.05, 
                         connectivity = NULL) {
    
    # checks
    if (!is.vector(data)) {
      stop("data needs to be a prepared vector excluding the 0 diagonal and including corresponding weights only once")
    }
    
    # resample from data
    resampled <- sample(data, size = n, replace = TRUE)
    
    # select only nonzeros
    nonzero <- which(resampled != 0)
    
    if (length(nonzero) > 0) {
      
      # always add jitter (cleaner logic)
      noise <- rnorm(length(nonzero), mean = 0, sd = jitter_sd)
      resampled[nonzero] <- resampled[nonzero] + noise
      
      resampled[nonzero] <- pmax(resampled[nonzero], min(data))
      resampled[nonzero] <- pmin(resampled[nonzero], max(data))
      
      # apply connectivity scaling if specified
      if (!is.null(connectivity)) {
        
        # check
        if (length(connectivity) != 1 || !connectivity %in% c("low", "moderate", "high")) {
          stop("connectivity should be one of: 'low', 'moderate', 'high'")
        }
        
        multiplier <- switch(connectivity,
                             "low" = 0.80,
                             "moderate" = 1.10,
                             "high" = 2)
        
        resampled[nonzero] <- resampled[nonzero] * multiplier
      }
    }
    
    return(resampled)
  }
  
  
  
  ## threshold function
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
  
  
  ## define base values 
  symptoms <- c("Depressed Mood", 
                "Pleasure/Interest",
                "Weight",
                "Insomnia/Hypersomnia", 
                "Psychomotor Agitation/Retardation",
                "Fatigue",
                "Worthlessness/Guilt",
                "Concentration",
                "Suicidality")
  
  index <- length(symptoms)
  
  low <- 0.08
  medium <- 1.10
  high <- 2.00
  
  
  # total number of weights needed
  n <- choose(index, 2) +     # weight number below the diagonal
    index                     # weight number on the diagonal
  
  
  
  ## generate weight matrix
  
  
  # calculate weights
  weight <- weight_sim(data = weight_data, 
                       n = n, 
                       jitter_sd = weight_jitter,
                       connectivity = connectivity)
  
  # create mirrored matrix filled with weights
  w <- matrix(0, index, index)
  w[upper.tri(w)] <- weight[1:choose(index, 2)]
  weights <- w + t(w)
  
  # if requested set diagonal (autocorrelations)
  if(auto == TRUE){
    diag(weights) <- weight[(choose(index, 2) + 1):
                              (choose(index, 2) + index)]
  }
  
  # assign symptom names
  rownames(weights) <- symptoms
  colnames(weights) <- symptoms
  
  
  ## generate threshold vector
  thresh <- thres_sim(data = thres_data,
                      n = index,
                      jitter_sd = thres_jitter)
  
  # assign symptom names
  names(thresh) <- symptoms
  
  
  
  ### create raw personalized parameter list
  params <- list(thresholds = thresh,
                 weights = weights)
  
  
  
  ##############################################################################
  ############################### LBA Maintenance ##############################
  
  # Depression:   depressive behavior; i.e., a general reduction in behavior
  
  
  # Maintenance Mechanisms
  
  # A:            negative reinforcement/lack of positive reinforcement lead to a 
  #               reduction in behavior which in turn reduces the probability of 
  #               experiencing positive reinforcement, hence, again reducing 
  #               behavior. 
  
  # B:            support first reinforces the depressed state due to increased 
  #               attention, and upon withdrawal leads to FEELING worse/desperate.
  #                 - although everything is defined purely BEHAVIORAL, here, 
  #                   FEELING - clearly not behavioral - is introduced. 
  
  
  # Network Setup
  
  # A:
  #   b:      (1) tendency to reduce behavior (in the absence of reinforcement)
  #               - reduced thresholds for behavioral (somatic) symptoms
  #           (2) responsiveness to external reinforcement initially remains
  
  #   w:      (3) both behavior and responsiveness to external reinforcement 
  #               influence each other strongly. behavior reduction leads to 
  #               reduced responsiveness and anhedonia first of all.
  #               - increased connectivity between behavioral symptoms and loss of
  #                 interest/pleasure
  
  # B:        (i) unclear as FEELING does not align with rest of model            ## double check this with Marie - I believe this is not suited and inappropriate
  #   b:      (1) potentially, tendency to feel guilty
  #               - reduced threshold for guilt/worthlessness
  
  #   w:      (2) potentially, guilt and behavioral symptoms influence each other 
  #               - increased connectivity between behavioral symptoms + guilt
  #
  #           (3) + connection guilt - suicidal
  #               - increased connectivity between suicidal and guilt
  
  ############################## b Formalization ###############################
  
  if(LBA_weights == TRUE){
                                                                                # should this be medium or high??????
    # adjust loss of pleasure - behavioral weights
    params$weights["Pleasure/Interest", "Fatigue"] <- params$weights["Pleasure/Interest", "Fatigue"] * high
    params$weights["Fatigue", "Pleasure/Interest"] <- params$weights["Pleasure/Interest", "Fatigue"] * high
    
    params$weights["Pleasure/Interest", "Weight"] <- params$weights["Pleasure/Interest", "Weight"] * high
    params$weights["Weight", "Pleasure/Interest"] <- params$weights["Pleasure/Interest", "Weight"] * high
    
    params$weights["Pleasure/Interest", "Psychomotor Agitation/Retardation"] <- params$weights["Pleasure/Interest", "Psychomotor Agitation/Retardation"] * high
    params$weights["Psychomotor Agitation/Retardation", "Pleasure/Interest"] <- params$weights["Pleasure/Interest", "Psychomotor Agitation/Retardation"] * high
    
    
    
    symptoms <- c("Depressed Mood", 
                  "Pleasure/Interest",
                  "Weight",
                  "Insomnia/Hypersomnia", 
                  "Psychomotor Agitation/Retardation",
                  "Fatigue",
                  "Worthlessness/Guilt",
                  "Concentration",
                  "Suicidality")
    
  }
  
  
  
  
  
  KBA_weights = FALSE,
  CT_weights = FALSE,
  
  thres_data,            # reference threshold generation data
  thres_jitter = 0.05,   # define how much thresholds should 
  
  # include CBT treatment maintenance thresholds
  LBA_thres = FALSE,
  KBA_thres = FALSE,
  CT_thres = FALSE,
  
  
  
  ### return dataframe as output
  return(pers.par)
}




################################################################################
################################ LBA Maintenance ###############################

# Depression:   depressive behavior; i.e., a general reduction in behavior


# Maintenance Mechanisms

# A:            negative reinforcement/lack of positive reinforcement lead to a 
#               reduction in behavior which in turn reduces the probability of 
#               experiencing positive reinforcement, hence, again reducing 
#               behavior. 

# B:            support first reinforces the depressed state due to increased 
#               attention, and upon withdrawal leads to FEELING worse/desperate.
#                 - although everything is defined purely BEHAVIORAL, here, 
#                   FEELING - clearly not behavioral - is introduced. 


# Network Setup

# A:
#   b:      (1) tendency to reduce behavior (in the absence of reinforcement)
#               - reduced thresholds for behavioral symptoms
#           (2) responsiveness to external reinforcement initially remains
         
#   w:      (3) both behavior and responsiveness to external reinforcement 
#               influence each other strongly. behavior reduction leads to 
#               reduced responsiveness and anhedonia first of all.
#               - increased connectivity between behavioral symptoms and loss of
#                 interest/pleasure

# B:        (i) unclear as FEELING does not align with rest of model            ## double check this with Marie - I believe this is not suited and inappropriate
#   b:      (1) potentially, tendency to feel guilty
#               - reduced threshold for guilt/worthlessness

#   w:      (2) potentially, guilt and behavioral symptoms influence each other 
#               - increased connectivity between behavioral symptoms + guilt
#
#           (3) + connection guilt - suicidal
#               - increased connectivity between suicidal and guilt



############################### b Formalization ################################
### function that adjust an individuals thresholds in line with network setup

LBA_thresh <- function(A = TRUE,     #  phenomena A will be implemented
                       B = FALSE){   #  phenomena B will be implemented
  
  low: c,l=0.80; medium: c,l=1.10; high: c,l=2.00; Cramer et al., 2016
  
}



################################################################################
################################ KBA Maintenance ###############################




################################################################################
################################ CT Maintenance ################################



################################################################################
########################### Treatment Formalisations ###########################
################################################################################     


         
         
         
         
         
         
         
         
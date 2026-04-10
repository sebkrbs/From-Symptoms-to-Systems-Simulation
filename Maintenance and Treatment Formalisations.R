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


# Network Phenomena

# A:
#   b:      (1) tendency to reduce behavior (in the absence of reinforcement)
#           (2) responsiveness to external reinforcement initially remains
         
#   w:      (3) both behavior and responsiveness to external reinforcement 
#               influence each other strongly. behavior reduction leads to 
#               reduced responsiveness and anhedonia first of all.

# B:        (i) unclear as FEELING does not align with rest of model            ## double check this with Marie - I believe this is not suited and inappropriate
#   b:      (1) potentially, tendency to feel guilty

#   w:      (2) potentially, guilt and behavioral symptoms influence each other 
#           (3) + connection guilt - suicidal.



############################### b Formalization ################################
LBA_thresh <- function(A = TRUE,     #  phenomena A will be implemented
                       B = FALSE){   #  phenomena B will be implemented
  
}



################################################################################
################################ KBA Maintenance ###############################




################################################################################
################################ CT Maintenance ################################



################################################################################
########################### Treatment Formalisations ###########################
################################################################################     


         
         
         
         
         
         
         
         
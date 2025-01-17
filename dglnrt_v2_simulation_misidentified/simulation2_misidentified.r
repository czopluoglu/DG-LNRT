require(MASS)
require(MBESS)
require(matrixStats)
require(cmdstanr)
require(here)
require(rstan)
require(psych)
###############################################################################

# This is a code being used for simulating item preknowledge data. 
# There are 1636 examinees and 170 items as in Real Dataset 2.
# There are 64 items compromised, and 46 examinees were assumed to have access
# to the 64 compromised items.

# DG-LNRT model is fitted as usual with an important difference. We assume that
# 50% of compromised items were incorrectly identified. So, while fitting the model
# we treated 64 items as compromised, but only 32 of these items were correctly
# identified as compromised and the remaining 32 items were incorrectly identified
# as compromised.
# In addition, 32 out of 106 uncompromised items were assumed to be incorrectly 
# identified as compromised, and they were treated as compromised items during 
# the model fitting process.


# Parameters come from EMIP paper where a multigroup lognormal response time
# model was fitted with gated mechanism using Real Dataset 2.

# For every replication, we generate a new set of model parameters from the 
# the defined distributions.

# Prepare an Appendix describing the model, model estimation, fitting, etc.
# Or, put the whole paper as an appendix. 

####################################################################
####################################################################
####################################################################

# A generic function to simulate DG-LNRT data with no item preknowledge

  sim_dglnrt <- function() {
    
    # MODEL PARAMETERS
    
    beta  <- rnorm(170,3.98,0.32)
    alpha <- rnorm(170,2.06,0.28)
    
    # Tau for unflagged examinees
    
    cor0 <- matrix(c(1,.95,.95,1),2,2)
    tau0 <- mvrnorm(1590,c(0,0),cor2cov(cor0,c(0.17,0.16)))
    
    # Tau for flagged examinees 
    
    cor1 <- matrix(c(1,.92,.92,1),2,2)
    tau1 <- mvrnorm(46,c(0.19,0.28),cor2cov(cor1,c(0.35,0.30)))
       
    # A vector for item status (0: not disclosed, 1:disclosed)
    
    C    <- c(rep(0,106),rep(1,64))
    
      # Shuffle the vector, randomly selected 64 items are compromised
    
        C <- C[base::sample(1:170,170)]
	
    # RESPONSE TIME GENERATION
    
    # Note that the response time data generated is
    # already on the log scale
    
    # Unflagged Examinees
    
    rt0 <- matrix(nrow = 1590, ncol = 170)
    
    for (i in 1:1590) {
      for (j in 1:170) {
        p_t = beta[j] - tau0[i, 1]
        p_c = beta[j] - tau0[i, 2]
        p   = p_t * (1 - C[j]) + p_c * C[j]
        rt0[i, j] = rnorm(1, p, 1 / alpha[j])
      }
    }
    
    # Flagged Examinees
    
    rt1 <- matrix(nrow = 46, ncol = 170)
    
    for (i in 1:46) {
      for (j in 1:170) {
        p_t = beta[j] - tau1[i, 1]
        p_c = beta[j] - tau1[i, 2]
        p   = p_t * (1 - C[j]) + p_c * C[j]
        rt1[i, j] = rnorm(1, p, 1 / alpha[j])
      }
    }
    
    # Combine the groups
    
    rt <- rbind(cbind(data.frame(exp(rt0)), gr = 1),
                cbind(data.frame(exp(rt1)), gr = 2))
    
    
    # Random shuffle of examinees
    
    rt <- rt[sample(1:nrow(rt),nrow(rt),replace = FALSE),]
    
    return(list(
      rt = rt,
      b = beta,
      a = alpha,
      tau_t = c(tau0[, 1], tau1[, 1]),
      tau_c = c(tau0[, 2], tau1[, 2]),
      C = C
    ))
  }
  
##############################################################################

# This section runs a single replication

  # Read the Stan model syntax, this is same across all replications
  
    mod <- cmdstan_model(here('dglnrt_v2_simulation_misidentified/dglnrt2.stan'))

  # Simulate data 
    
    data <- sim_dglnrt()
    
    d.sub <- data$rt[,1:170]
    
    d.sub$ID <- 1:1636
   
    
  # Rehape the data into the long format
    
    d.long <- reshape(
      data        = d.sub,
      idvar       = "ID",
      varying     = list(colnames(d.sub)[1:170]),
      timevar     = "Item",
      times       = 1:170,
      v.names      = "RT",
      direction   = "long"
    )
    
    d.long <- na.omit(d.long)
    d.long$logRT <- log(d.long$RT)
    
    d.long$i.status <- NA
    
    # Randomly select 32 out of 64 compromised items
    # and feed them as compromised items
    # the remaining 32 compromised items will be assumed not detected, 
    # treated as normal during the model fitting process
    
    # In addition, randomly select another 32 items from uncompromised items
    # and treat them as compromised during the model fitting process
    
    
    set1 <- base::sample(which(data$C==1),32) 
    set2 <- base::sample(which(data$C==0),32) 
    
    flagged.item <- c(set1,set2)
    
    
    d.long[which(d.long$Item%in%flagged.item==TRUE),]$i.status = 1
    d.long[which(d.long$Item%in%flagged.item==FALSE),]$i.status = 0
    
    
  # Data object for Stan
    
    data_rt <- list(
      J              = 170,
      I              = 1636,
      n_obs          = length(d.long$RT),
      ind_person_obs = d.long$ID,
      ind_item_obs   = d.long$Item,
      i_status_obs   = d.long$i.status,
      Y              = d.long$logRT
    )
    
    
  # Fit the model using cmdstan
    
    fit <- mod$sample(
      data = data_rt,
      seed = 1234,
      chains = 4,
      parallel_chains = 4,
      iter_warmup   = 1000,
      iter_sampling = 2000,
      refresh = 100,
      adapt_delta = 0.99)

    
  # Save the output file   
    
    fit$cmdstan_summary()
    
    stanfit <- rstan::read_stan_csv(fit$output_files())
    
    
    save.image(here(''))

    
    # 100 replications were run independently in the UO computing cluster (Talapas)
    # due to the computational demand 
    
    # See the code under /dglnrt_v2_simulation/talapas for more syntax and specific 
    # files

################################################################################
################################################################################
################################################################################
#
# Outcome Analysis
#
################################################################################
################################################################################
################################################################################
    
# Read the output files for 100 replications
    
# Create two list objects with length of 100 to save the output for each replication
    
f <- list.files(here('data/dglnrt_v2_simulation_misidentified'))

stanfit.list <- vector('list',100)
data.list    <- vector('list',100)

for(kk in 1:100){
  
  ch <- file.exists(here(paste0('data/dglnrt_v2_simulation_misidentified/rep',kk,'.RData')))  
  
  if(ch==TRUE){
    
    load(here(paste0('data/dglnrt_v2_simulation_misidentified/rep',kk,'.RData')))
    stanfit.list[[kk]] <- stanfit
    data.list[[kk]]    <- data
    print(kk)
    
    rm(list = ls()[!ls()%in%c('stanfit.list','data.list')])
  }

}

################################################################################
# For each replication, extract the estimate of T for each individual
# Save them in a 100 x 1636 matrix
# Each row represents a replication
# Each column represents an individual within a replication
# Cell values are the estimate of posterior probability of item preknowledge
# for an individual in a replication


param <- matrix(nrow=100,ncol=1636)

for(i in 1:100){
  
  if(is.null(stanfit.list[[i]])==FALSE){
    
    fit   <- stanfit.list[[i]]
    Ts <- as.numeric(summary(fit, pars = c("T"), probs = c(0.025, 0.975))$summary[,1])
    param[i,] = Ts
  }
  
  print(i)
}

param <- param[1:97,]


# For a given cut-off value, compute the average proportion of falsely 
# identified individuals across 100 replications


th = 0.999

out <- data.frame(matrix(NA,100,4))
colnames(out) <- c('FP','TP')

for(R in 1:100){
  
  Ts   <- param[R,]
  t    <- ifelse(Ts>th,1,0)
  true <- ifelse(data.list[[R]]$rt$gr==2,1,0)
  tab  <- table(t,true)
  
  out[R,]$FP  = sum(t==1 & true==0)
  out[R,]$TP  = sum(t==1 & true==1)
}

# FPR across 100 replications

sum(out$FP)/(1590*100)

# TPR across 100 replications

sum(out$TP)/(46*100)

# Precision across 100 replications

sum(out$TP)/(sum(out$FP)+sum(out$TP))

################################################################################
# Check item parameter recovery across 100 replications
################################################################################

# Betas


tr <- matrix(NA,100,170)
est <- matrix(NA,100,170)

for(i in 1:97){
  fit     <- stanfit.list[[i]]
  est[i,] <- as.numeric(summary(fit, pars = c("beta"), probs = c(0.025, 0.975))$summary[,1])
  tr[i,]  <-  data.list[[i]]$b
  print(i)
}

tr  <- tr[1:97,]
est <- est[1:97,] 

bias1 <- c()
bias2 <- c()

for(i in 1:97){
  
  temp <- tr[i,] - est[i,]
  
  bias1[i] = mean(temp[which(data.list[[i]]$C==0)])
  
  bias2[i] = mean(temp[which(data.list[[i]]$C==1)])
  
}

round(mean(bias1),3)
round(mean(bias2),3)

# Alphas

tr <- matrix(NA,100,170)
est <- matrix(NA,100,170)

for(i in 1:97){
  fit     <- stanfit.list[[i]]
  est[i,] <- as.numeric(summary(fit, pars = c("alpha"), probs = c(0.025, 0.975))$summary[,1])
  tr[i,]  <-  data.list[[i]]$a
  print(i)
}

tr  <- tr[1:97,]
est <- est[1:97,] 

bias1 <- c()
bias2 <- c()

for(i in 1:97){
  
  temp <- tr[i,] - est[i,]
  
  bias1[i] = mean(temp[which(data.list[[i]]$C==0)])
  
  bias2[i] = mean(temp[which(data.list[[i]]$C==1)])
  
}

round(mean(bias1),3)
round(mean(bias2),3)











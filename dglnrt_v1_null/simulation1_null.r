require(MASS)
require(MBESS)
require(matrixStats)
require(cmdstanr)
require(here)
require(rstan)
require(psych)
###############################################################################

# This is a code being used for simulating a null data (where there is no item 
# preknowledge). There are 93 examinees and 25 items as in Real Dataset 1.
# DG-LNRT model is fitted as usual. We expect DG-LNRT to detect nobody in this 
# dataset except false positives.

# Parameters come from EMIP paper where a multigroup lognormal response time
# model was fitted with gated mechanism using Real Dataset 1.

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
    
    beta  <- rnorm(25,4.19,0.38)
    alpha <- rnorm(25,1.42,0.37)
    
    # Tau for control group
    
    cor0 <- matrix(c(1,.84,.84,1),2,2)
    tau0 <- mvrnorm(33,c(0,0),cor2cov(cor0,c(0.18,0.17)))
    
    # Tau for experimental group (Items Disclosed)
    
    cor1 <- matrix(c(1,.56,.56,1),2,2)
    tau1 <- mvrnorm(30,c(0.05,1.34),cor2cov(cor1,c(0.55,0.47)))
    
    # Tau for experimental group (Items and Answers Disclosed)
    
    cor2 <- matrix(c(1,.43,.43,1),2,2)
    tau2 <- mvrnorm(30,c(-0.18,1.47),cor2cov(cor2,c(0.43,0.46)))
    
    # A vector for item status (0: not disclosed, 1:disclosed)
    
    C    <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
    
    # RESPONSE TIME GENERATION
    
    # Note that the response time data generated is
    # already on the log scale
    
    # Control Group
    
    rt0 <- matrix(nrow = 33, ncol = 25)
    
    for (i in 1:33) {
      for (j in 1:25) {
        p_t = beta[j] - tau0[i, 1]
        p_c = beta[j] - tau0[i, 2]
        p   = p_t * (1 - C[j]) + p_c * C[j]
        rt0[i, j] = rnorm(1, p, 1 / alpha[j])
      }
    }
    
    # Experimental Group 1 (Item Disclosed)
    
    rt1 <- matrix(nrow = 30, ncol = 25)
    
    for (i in 1:30) {
      for (j in 1:25) {
        p_t = beta[j] - tau1[i, 1]
        p_c = beta[j] - tau1[i, 2]
        p   = p_t * (1 - C[j]) + p_c * C[j]
        rt1[i, j] = rnorm(1, p, 1 / alpha[j])
      }
    }
    
    # Experimental Group 2 (Item and Answers Disclosed)
    
    rt2 <- matrix(nrow = 30, ncol = 25)
    
    for (i in 1:30) {
      for (j in 1:25) {
        p_t = beta[j] - tau2[i, 1]
        p_c = beta[j] - tau2[i, 2]
        p   = p_t * (1 - C[j]) + p_c * C[j]
        rt2[i, j] = rnorm(1, p, 1 / alpha[j])
      }
    }
    
    # Combine the groups
    
    rt <- rbind(cbind(data.frame(exp(rt0)), gr = 1),
                cbind(data.frame(exp(rt1)), gr = 2),
                cbind(data.frame(exp(rt2)), gr = 3))
    
    
    return(list(
      rt = rt,
      b = beta,
      a = alpha,
      tau_t = c(tau0[, 1], tau1[, 1], tau2[,1]),
      tau_c = c(tau0[, 2], tau1[, 2], tau2[,2])
    ))
  }

##############################################################################

# RUN 100 Replications
# Simulate data and fit DG-LNRT

  sim.datasets <- vector('list',100)
  stanfit <- vector('list',100)
  
  mod <- cmdstan_model(here('dglnrt_v1_null/dglnrt2.stan'))
  
  for(i in 1:100){
    
    data <- sim_dglnrt()
    
    sim.datasets[[i]] = data
    
    d.sub <- data$rt[,1:25]
    
    d.sub$ID <- 1:93
    
    d.long <- reshape(
      data        = d.sub,
      idvar       = "ID",
      varying     = list(colnames(d.sub)[1:25]),
      timevar     = "Item",
      times       = 1:25,
      v.names      = "RT",
      direction   = "long"
    )
    
    d.long$i.status <- ifelse(d.long$Item%%2==0,1,0)
    
    d.long$logRT <- log(d.long$RT)
    
    data_rt <- list(
      J              = 25,
      I              = 93,
      n_obs          = length(d.long$RT),
      ind_person_obs = d.long$ID,
      ind_item_obs   = d.long$Item,
      i_status_obs   = d.long$i.status,
      Y              = d.long$logRT
    )
    
    
    fit <- mod$sample(
      data = data_rt,
      seed = 1234,
      chains = 4,
      parallel_chains = 4,
      iter_warmup   = 1000,
      iter_sampling = 2000,
      refresh = 1000,
      adapt_delta = 0.99
    )
    
    fit$cmdstan_summary()
    
    stanfit[[i]] <- rstan::read_stan_csv(fit$output_files())
  }  
###############################################################################


save.image(here('dglnrt_v1_null/dglnrt_v1_null.RData'))

###################################################################################

################################################################################
################################################################################
################################################################################
#
# Outcome Analysis
#
################################################################################
################################################################################
################################################################################
################################################################################
  
  load(here("data/dglnrt_v1_null/dglnrt_v1_null.RData"))
  
  
# T
  
  param <- matrix(nrow=100,ncol=93)
  corr <- c()
  
  for(i in 1:100){
    fit   <- stanfit[[i]]
    Ts <- as.numeric(summary(fit, pars = c("T"), probs = c(0.025, 0.975))$summary[,1])
    param[i,] = Ts
  }
  
  hist(param)
  
  th = 0.99
  
  round(c(mean(rowMeans(param>th)),
          min(rowMeans(param>th)),
          max(rowMeans(param>th))),3)
  
  
# Betas
  
  param <- matrix(nrow=100,ncol=25)
  corr <- c()
  
  for(i in 1:100){
    fit   <- stanfit[[i]]
    betas <- as.numeric(summary(fit, pars = c("beta"), probs = c(0.025, 0.975))$summary[,1])
    param[i,] = betas - sim.datasets[[i]]$b
    corr[i] = cor(betas, sim.datasets[[i]]$b)
    print(i)
  }
  
  mean(rowMeans(param)) # average bias
  min(rowMeans(param)) # average bias
  max(rowMeans(param)) # average bias
  
  mean(sqrt(rowMeans(param^2))) # average root mean squared error
  min(sqrt(rowMeans(param^2))) # average root mean squared error
  max(sqrt(rowMeans(param^2))) # average root mean squared error
  
   
  mean(corr)            # average correlation
  min(corr)
  max(corr)
  
# Alphas
  
  param <- matrix(nrow=100,ncol=25)
  corr <- c()
  
  for(i in 1:100){
    fit   <- stanfit[[i]]
    alphas <- as.numeric(summary(fit, pars = c("alpha"), probs = c(0.025, 0.975))$summary[,1])
    param[i,] = alphas - sim.datasets[[i]]$a
    corr[i] = cor(alphas, sim.datasets[[i]]$a)
    print(i)
  }
  

  mean(rowMeans(param)) # average bias
  min(rowMeans(param)) # average bias
  max(rowMeans(param)) # average bias
  
  mean(sqrt(rowMeans(param^2))) # average root mean squared error
  min(sqrt(rowMeans(param^2))) # average root mean squared error
  max(sqrt(rowMeans(param^2))) # average root mean squared error
  
  
  mean(corr)            # average correlation
  min(corr)
  max(corr)
  
# Tau_t
  
  param <- matrix(nrow=100,ncol=93)
  corr <- c()
  
  for(i in 1:100){
    fit   <- stanfit[[i]]
    taus <- as.numeric(summary(fit, pars = c("tau_t"), probs = c(0.025, 0.975))$summary[,1])
    param[i,] = taus - sim.datasets[[i]]$tau_t
    corr[i] = cor(taus, sim.datasets[[i]]$tau_t)
  }
  
  mean(rowMeans(param)) # average bias
  
  sqrt(mean(param^2))   # root mean squared error
  
  mean(corr)            # average correlation

# Tau_c
  
  param <- matrix(nrow=100,ncol=93)
  corr <- c()
  
  for(i in 1:100){
    fit   <- stanfit[[i]]
    taucs <- as.numeric(summary(fit, pars = c("tau_c"), probs = c(0.025, 0.975))$summary[,1])
    param[i,] = taucs - sim.datasets[[i]]$tau_c
    corr[i] = cor(taucs, sim.datasets[[i]]$tau_c)
  }
  
  mean(rowMeans(param)) # average bias
  
  sqrt(mean(param^2))   # root mean squared error
  
  mean(corr)            # average correlation









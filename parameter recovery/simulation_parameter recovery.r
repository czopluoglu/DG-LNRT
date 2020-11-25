

require(MASS)
require(MBESS)
require(matrixStats)
require(cmdstanr)
require(here)
require(rstan)

###############################################################################

# This is the code being used for a simulation study to examine the quality of 
# item parameter estimates relevant to Real Dataset 1 when sample size is 93 and
# number of items is 25, and when the responses for 60 examinees are missing for 
# the same set of 12 items. 

# This is mimicking the Real Dataset 1 with no contamination of item preknowledge). 
# The distributions of $\beta$, $\alpha$, and $\tau$ were slightly different 
# than the ones observed in the Real Dataset 1.

#############################################################################

# A generic function to simulate LNRT data.


sim_dglnrt <- function() {

  # MODEL PARAMETERS

    # For each replication, the beta, alpha, and tau were re-generated based on the
    # following distributions
    
    # beta ~ N(4.5,0.5)
    # alpha ~ N(1.2,0.4)
    # tau ~ 0, .5
  
    beta  <- rnorm(25,4.5,0.5)
    alpha <- (rnorm(25,1.2,0.4))
      while(sum(alpha<.3) > 0){
        alpha <- (rnorm(25,1.2,0.4))
      }

    tau <- rnorm(93,0,.5)
    
  # RESPONSE TIME DATA GENERATION

    rt <- matrix(nrow = length(tau), ncol = length(beta))
    
    for (i in 1:length(tau)) {
      for (j in 1:length(beta)) {
        rt[i, j] = rnorm(1,beta[j]-tau[i], 1 / alpha[j])
      }
    }
    
    return(list(
      rt = rt,
      b = beta,
      a = alpha,
      t = tau
    ))
    
}

##############################################################################

# List objects to store the datasets and fitted model objects

  sim.datasets <- vector('list',100)
  stanfit <- vector('list',100)

# Stan model syntax
  
  mod <- cmdstan_model(here('parameter recovery/dglnrt_v1_parameter_recovery.stan'))

# Run 100 replications
    
    for(i in 1:100){
      
      data <- sim_dglnrt()
      
      sim.datasets[[i]] = data
      
      d.sub <- as.data.frame(data$rt)
    
      d.sub$ID <- 1:93
      
      # Assign missing value for the responses of 60 examinees on 12 items
      
      d.sub[34:93,seq(2,25,2)]=NA
    
      # Tranform to a long format data
      
      d.long <- reshape(
        data        = d.sub,
        idvar       = "ID",
        varying     = list(colnames(d.sub)[1:25]),
        timevar     = "Item",
        times       = 1:25,
        v.names      = "RT",
        direction   = "long"
      )
      
      # remove missing data 
      
      d.long <- na.omit(d.long)
      
      # data list object for stan
      
      data_rt <- list(
        J              = 25,
        I              = 93,
        n_obs          = length(d.long$RT),
        ind_person_obs = d.long$ID,
        ind_item_obs   = d.long$Item,
        Y              = d.long$RT
      )
      
      # fit the model
      
      fit <- mod$sample(
        data = data_rt,
        seed = 1234,
        chains = 4,
        parallel_chains = 4,
        iter_warmup   = 1000,
        iter_sampling = 2000,
        refresh = 100,
        adapt_delta = 0.99
      )
      
      # Extract the output files
    
      fit$cmdstan_summary()
      
      stanfit[[i]] <- rstan::read_stan_csv(fit$output_files())
    }  
###############################################################################

# Store everything

  save.image(here('parameter recovery/parameter recovery_N=93.RData'))

###################################################################################

# OUTCOME ANALYSiS  
  
param <- matrix(nrow=100,ncol=55)


for(i in 1:100){
  
  fit <- stanfit[[i]]

  param[i,1:3] = as.numeric(summary(fit, pars = c("mu1","sigma1","sigma_t"), probs = c(0.025, 0.975))$summary[,1])
  
  alphas <- summary(fit, pars = c("alpha"), probs = c(0.025, 0.975))$summary[,1]
  param[i,4:5]   = c(mean(alphas),sd(alphas))
  print(i)
}

param <- matrix(nrow=100,ncol=25)
corr <- c()

for(i in 1:100){
  fit   <- stanfit[[i]]
  betas <- as.numeric(summary(fit, pars = c("beta"), probs = c(0.025, 0.975))$summary[,1])
  param[i,] = betas - sim.datasets[[i]]$b
  corr[i] = cor(betas, sim.datasets[[i]]$b)
}

colMeans(param)
sqrt(colSums(param^2)/100)
mean(corr)


param <- matrix(nrow=100,ncol=25)
corr <- c()

for(i in 1:100){
  fit   <- stanfit[[i]]
  alphas <- as.numeric(summary(fit, pars = c("alpha"), probs = c(0.025, 0.975))$summary[,1])
  param[i,] = alphas - sim.datasets[[i]]$a
  corr[i] = cor(alphas, sim.datasets[[i]]$a)
}

colMeans(param)
sqrt(colSums(param^2)/100)
mean(corr)


param <- matrix(nrow=100,ncol=93)
corr <- c()

for(i in 1:100){
  fit   <- stanfit[[i]]
  taus <- as.numeric(summary(fit, pars = c("tau_t"), probs = c(0.025, 0.975))$summary[,1])
  param[i,] = taus - sim.datasets[[i]]$t
  corr[i] = cor(taus, sim.datasets[[i]]$t)
}

colMeans(param)
sqrt(colSums(param^2)/100)
mean(corr)










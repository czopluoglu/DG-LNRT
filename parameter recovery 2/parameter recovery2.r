
require(MASS)
require(MBESS)
require(matrixStats)
require(cmdstanr)
require(here)
require(rstan)



################################################################################
################################################################################
################################################################################

# This folder includes the code being used for a simulation study to examine the 
# quality of item parameter estimates relevant to Real Dataset 1  when sample 
# size is 93 and number of items is 25, and when the responses for 60 examinees 
# are missing for the same 12 items. 

# First, the LNRT model was fitted to uncontaminated data after excluding the
# responses of 60 examinees to 12 compromised items. Then, the estimated 
# model parameters were saved.

# Second, a simulation study is run by mimicking the Real Dataset 1 setting 
# with no contamination of item preknowledge). The exact same set of $\beta$, 
# $\alpha$, and $\tau$ parameters estimated from data were being used and 
# fixed across 100 replications.

################################################################################
################################################################################
################################################################################
#
# Fitting model to uncontaminated Real Dataset 1 
#
################################################################################
################################################################################
################################################################################
################################################################################

# Import data

  d <- read.csv(here('data/uva_rt.csv'),na.strings = '#N/A')


  d.sub <- d[,c("ï..ID","COND","Q1RT","Q2RT","Q3RT","Q4RT","Q5RT","Q6RT","Q7RT",
                "Q8RT","Q9RT","Q10RT","Q11RT","Q12RT","Q13RT","Q14RT","Q15RT",
                "Q16RT","Q17RT","Q18RT","Q19RT","Q20RT","Q21RT","Q22RT","Q23RT",
                "Q24RT","Q25RT")]

  colnames(d.sub)[1] = 'ID'

  d.sub$ID <- 1:93

# Remove the responses for 12 disclosed items 
# from individuals in Experimental Group 1 and 2
# This leaves a big chunk of missingness in the dataset

  d.sub[which(d.sub$COND==2 | d.sub$COND==3),seq(4,27,2)] = NA

# Transform to long format
  
  d.long <- reshape(data        = d.sub,
                    idvar       = "ID",
                    varying     = list(colnames(d.sub)[3:27]),
                    timevar     = "Item",
                    times       = 1:25,
                    v.names      = "RT",
                    direction   = "long")
# Remove missing data
  
  d.long <- na.omit(d.long)

# Take log of response times
  
  d.long$logRT <- log(d.long$RT)

# Data list for Stan 
  
  data_rt <- list(J              = 25,
                  I              = 93,
                  n_obs          = length(d.long$RT),
                  ind_person_obs = d.long$ID,
                  ind_item_obs   = d.long$Item,
                  Y              = d.long$logRT)

##########################################################

# Below is the logic and references used for PRIORS

# van der Linder 2006
# page 191, Eq 23-24

# Best estimates of average time intensity and time-discrimination

  beta0 = colMeans(log(d.sub[,3:27]),na.rm=TRUE)

  alpha0 = sqrt(1/colMeans((log(d.sub[,3:27]) - matrix(beta0,93,25,byrow=T))^2,na.rm=TRUE))

  mean(beta0)
  mean(alpha0)

  N= 93
  a = mean(alpha0)

  # hyper parameters for inverse gamma 
  # Levy, Bayesian Psychometric Analysis, page 83, inverse gamma prior distribution

  v0 = 93/2
  v0
  v0*a

  # alpha ~ InvGamm(46.5, 71.96)

require(invgamma)

test.gamma = rinvgamma(10000,46.5,71.96)
hist(test.gamma)
mean(test.gamma)

test.gamma = rgamma(10000,46.5,71.96)
hist(test.gamma)
mean(test.gamma)


# Statistical Rethinking

# Prior choice for sigma

# The exponential prior—dexp(1) in R code—has a much thinner tail than the Cauchy does. This
# induces more conservatism in estimates and can help your Markov chain converge correctly. The
# exponential is also the maximum entropy prior for the standard deviation, provided all we want to
# say a priori is the expected value. That is to say that the only information contained in an exponential
# prior is the mean value and the positive constraint.

# page 452, Gelman (2006) recommends the half-Cauchy because it is approximately uniform in the tail and 
# still weak near zero, without the odd behavior of traditional priors like the inverse-gamma. 
# Polson and Scott (2012) examine this prior in more detail. Simpson et al. (2014) also note that 
# the half-Cauchy prior has useful features, but recommend instead an exponential prior.
# Either is equally useful in all the examples in this book. [249]



##########################################################

# Stan Model Syntax

  mod <- cmdstan_model(here('parameter_recovery_2/dglnrt2.stan'))

# Fit the model

  fit <- mod$sample(
    data = data_rt,
    seed = 1234,
    chains = 4,
    parallel_chains = 4,
    iter_warmup   = 5000,
    iter_sampling = 10000,
    refresh = 100
  )

# Extract the output files
  
  fit$cmdstan_summary()

  stanfit <- rstan::read_stan_csv(fit$output_files())

# Save the information from fitting the model to uncontaminated Real Dataset 1
  
  save.image(here('data/simulation1/original results.RData'))

# Extract the model parameter estimates
  
  summary(stanfit, pars = c("mu1","sigma1","sigma_t"), 
          probs = c(0.025, 0.975))$summary

  betas <- summary(stanfit, pars = c("beta"), probs = c(0.025, 0.975))$summary
  betas 
  describe(betas[,1])

  alphas <- summary(stanfit, pars = c("alpha"), probs = c(0.025, 0.975))$summary
  alphas
  describe(alphas[,1])

  tau_t <- summary(stanfit, pars = c("tau_t"), probs = c(0.025, 0.975))$summary
  tau_t
  describe(tau_t[,1])

################################################################################
################################################################################
################################################################################
#
# Fitting model to uncontaminated Real Dataset 1 
#
################################################################################
################################################################################
################################################################################
################################################################################

# True model parameters to be used for simulation.
# These parameters are fixed across replications
  
true.beta <- c(4.068,3.999,4.695,3.615,4.659,4.032,4.268,4.112,4.848,4.360,
               5.133,4.071,4.381,4.778,4.577,4.694,3.730,3.771,4.455,3.644,
               3.818,4.069,3.893,4.046,4.201)

true.alpha <- c(1.796,1.974,1.919,1.565,1.874,2.033,1.699,1.627,1.758,1.757,
                2.090,1.661,1.516,1.250,1.358,1.863,1.578,2.039,1.353,1.965,
                1.195,1.580,1.063,1.032,0.945)

true.tau   <- c(-0.044,-0.274,-0.424,-0.416,-0.292,0.191,-0.005,0.032,-0.014,
                0.227,0.297,-0.022,-0.178,-0.310,-0.075,-0.107,0.034,0.833,
                -0.286,-0.119,-0.115,0.076,-0.138,-0.035,-0.187,-0.138,0.021,
                0.127,-0.289,0.099,-0.231,1.509,-0.114,-0.194,0.284,0.183,0.154,
                -0.061,-0.051,0.033,-0.162,0.087,-0.148,-0.369,0.030,-0.172,
                -0.426,0.848,-0.127,0.261,-0.395,0.194,0.034,-0.050,-0.445,
                -0.459,-0.001,1.818,0.177,0.075,-0.053,-0.360,-0.245,0.266,
                0.134,-0.094,0.073,0.646,0.125,-0.453,0.362,-0.008,-0.172,0.289,
                -0.253,-0.196,0.108,0.147,-0.187,0.015,0.738,-0.001,-0.478,0.154,
                0.225,0.012,-0.009,-0.422,-0.129,-0.242,-0.213,-0.445,-0.096)



################################################################################
################################################################################
################################################################################
#
# Simulation Study 
#
################################################################################
################################################################################
################################################################################
################################################################################

# A generic function to simulate LNRT data given the true model parameters

  sim_dglnrt <- function(beta, alpha, tau) {
    rt <- matrix(nrow = length(tau), ncol = length(beta))
    
    for (i in 1:length(tau)) {
      for (j in 1:length(beta)) {
        rt[i, j] = rnorm(1, beta[j] - tau[i], 1 / alpha[j])
      }
    }
  
    return(list(
      rt = rt,
      b = beta,
      a = alpha,
      t = tau
    ))
  }

# List objects to store the simulated datasets and fitted objects

  sim.datasets <- vector('list',100)
  stanfit <- vector('list',100)

# Stan Model Syntax
  mod <- cmdstan_model(here('parameter_recovery_2/dglnrt2.stan'))

# Run 100 replications, generate data and then fit the model
    
  for(i in 1:100){
    
    data <- sim_dglnrt(beta = true.beta, alpha = true.alpha, tau=true.tau)
    
    sim.datasets[[i]] = data
    
    d.sub <- as.data.frame(data$rt)
    d.sub[34:93,seq(2,25,2)]=NA
    
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
   
    d.long <- na.omit(d.long)
    
   
    data_rt <- list(
      J              = 25,
      I              = 93,
      n_obs          = length(d.long$RT),
      ind_person_obs = d.long$ID,
      ind_item_obs   = d.long$Item,
      Y              = d.long$RT
    )
    
    
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
  
    fit$cmdstan_summary()
    
    stanfit[[i]] <- rstan::read_stan_csv(fit$output_files())
  }  

# Store all information 
  
  save.image(here('parameter_recovery_2/parameter_recovery.RData'))

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
  




















################################################################################

# This fits the model to the real dataset by Toton and Maynes (2019)

################################################################################

require(rstan)
require(cmdstanr)
require(here)
require(psych)

################################################################################

# Import the dataset from Sarah and Toton 2019

d <- read.csv(here('data/uva_data_RA.csv'),na.strings = '#N/A')

# Update IDs such that ID goes from 1 to 93

colnames(d.sub)[1] <- 'ID'

d$ID <- 1:93

# Reshape data from wide format to long format

d.long <- reshape(data        = d,
                  idvar       = "ID",
                  varying     = list(colnames(d)[3:27]),
                  timevar     = "Item",
                  times       = 1:25,
                  v.names      = "score",
                  direction   = "long")

# Create a variable for the item compromised status. 
# Even-numbered items were disclosed and odd-numbered items were not disclosed

d.long$i.status <- ifelse(d.long$Item%%2==0,1,0)

# Remove missing data in the long format.
# If someone didn't respond to an item, only that particular item was removed
# the rest of the data for other items for the person remained.

d.long <- na.omit(d.long)

# data object for the Stan model above.

data_ra <- list(J              = 25,
                I              = 93,
                n_obs          = length(d.long$score),
                ind_person_obs = d.long$ID,
                ind_item_obs   = d.long$Item,
                i_status_obs   = d.long$i.status,
                Y              = d.long$score)
################################################################################
# STAN model syntax for the DG-IRT model
    
'    
data{
  int <lower=1> I;                           // number of examinees
  int <lower=1> J;                           // number of items
  int <lower=1> n_obs;                       // number of observations
  int <lower=1> ind_person_obs[n_obs];       // examinee position index
  int <lower=1> ind_item_obs[n_obs];         // item position index
  int <lower=0,upper=1> i_status_obs[n_obs]; // item status index
  int <lower=0,upper=1> Y[n_obs];            // observed binary responses
}

parameters {
 real mu_thetat;             // mean for theta_t
 real<lower=0> sigma_thetat; // sd for theta_t
  
 real mu_thetac;             // mean for theta_c
 real<lower=0> sigma_thetac; // sd for theta_c

 vector[I] theta_t;           // vector of true theta parameters
 vector[I] theta_c;           // vector of cheating theta parameters
 vector[J] b;                // vector of item difficulty parameters

}


transformed parameters {
  vector[I] T;
  
  for (i in 1:n_obs) {
    if(theta_t[ind_person_obs[i]]>theta_c[ind_person_obs[i]])
      
      T[ind_person_obs[i]] = 0;
    
    else 
      
      T[ind_person_obs[i]] = 1;
  }
}


model{
  
  sigma_thetat ~ exponential(1);
  sigma_thetac ~ exponential(1);
  mu_thetat    ~ normal(0,10);
  mu_thetac    ~ normal(0,10);

  theta_t    ~ normal(mu_thetat,sigma_thetat);
  theta_c    ~ normal(mu_thetac,sigma_thetac);

  b          ~ normal(0,1);

  for (i in 1:n_obs) {
    
    real p_t = inv_logit(theta_t[ind_person_obs[i]] - b[ind_item_obs[i]]);
    real p_c = inv_logit(theta_c[ind_person_obs[i]] - b[ind_item_obs[i]]);
    
    real p = (p_t^(1-T[ind_person_obs[i]]))*
      (((1-i_status_obs[i])*p_t + 
          (i_status_obs[i])*p_c)^T[ind_person_obs[i]]);
    
    Y[i] ~ bernoulli(p);
  }
}
'
    
###############################################################################

# Read the Stan Model Syntax 
    
mod <- cmdstan_model(here('dgirt_v1/dgirt.stan'))

# Fit the model using cmdstan
    
fit <- mod$sample(data = data_ra,
                  seed = 1234,
                  chains = 4,
                  parallel_chains = 4,
                  iter_warmup   = 1000,
                  iter_sampling = 4000,
                  refresh = 10
                  )

# Save the output object

fit$cmdstan_summary()

stanfit <- rstan::read_stan_csv(fit$output_files())

save.image(here('data/dgirt_v1/results_dgirt_v1.RData'))

################################################################################
#                               OUTPUT ANALYSIS

# Load the pre-saved model output

load(here('data/dgirt_v1/results_dgirt_v1.RData'))

# Summary of all parameters

params <- summary(stanfit, 
                  pars = c("mu_thetat","sigma_thetat",
                           "mu_thetac","sigma_thetac"),
                  probs = c(0.025, 0.975))$summary

params
# Extract the b parameters

b <- summary(stanfit, pars = c("b"), probs = c(0.025, 0.975))$summary
b
describe(b[,1])

# Extract the theta_t

theta_t <- summary(stanfit, pars = c("theta_t"), probs = c(0.025, 0.975))$summary
theta_t
describe(theta_t[,1])

# Extract the theta_c

theta_c <- summary(stanfit, pars = c("theta_c"), probs = c(0.025, 0.975))$summary
theta_c
describe(theta_c[,1])

# Extract Ts

Ts <- summary(stanfit, pars = c("T"), probs = c(0.025, 0.975))$summary
Ts

# For a given threshold, check the number of identified individuals in the data 
# and compare it to their known status of item preknowledge

th = .99

table(ifelse(Ts[,1]>th,1,0),d$COND)

 # FP: 0/33
 # TP: 7/60
 # PR: 7/7


hist(Ts[,1])


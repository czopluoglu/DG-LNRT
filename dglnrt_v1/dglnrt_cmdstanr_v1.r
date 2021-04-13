
################################################################################

# This fits the model to the real dataset by Toton and Maynes (2019)

################################################################################

require(cmdstanr)
require(here)
require(psych)

################################################################################

# Import the dataset from Sarah and Toton 2019

d <- read.csv(here('data/uva_rt.csv'),na.strings = '#N/A')


# Create a subset by select the response time data

d.sub <- d[,c("ï..ID","COND","Q1RT","Q2RT","Q3RT","Q4RT","Q5RT","Q6RT","Q7RT",
              "Q8RT","Q9RT","Q10RT","Q11RT","Q12RT","Q13RT","Q14RT","Q15RT",
              "Q16RT","Q17RT","Q18RT","Q19RT","Q20RT","Q21RT","Q22RT","Q23RT",
              "Q24RT","Q25RT")]

# Update IDs such that ID goes from 1 to 93

colnames(d.sub)[1] <- 'ID'

d.sub$ID <- 1:93

# Reshape data from wide format to long format

d.long <- reshape(data        = d.sub,
                  idvar       = "ID",
                  varying     = list(colnames(d.sub)[3:27]),
                  timevar     = "Item",
                  times       = 1:25,
                  v.names      = "RT",
                  direction   = "long")

# Create a variable for the item compromised status. 
# Even-numbered items were disclosed and odd-numbered items were not disclosed

d.long$i.status <- ifelse(d.long$Item%%2==0,1,0)

# Remove missing data in the long format.
# If someone didn't respond to an item, only that particular item was removed
# the rest of the data for other items for the person remained.

d.long <- na.omit(d.long)

# Transform the response time to log scale

d.long$logRT <- log(d.long$RT)

# data object for the Stan model above.

data_rt <- list(J              = 25,
                I              = 93,
                n_obs          = length(d.long$RT),
                ind_person_obs = d.long$ID,
                ind_item_obs   = d.long$Item,
                i_status_obs   = d.long$i.status,
                Y              = d.long$logRT)

##########################################################

# PRIORS


# van der Linder 2006
# page 191, Eq 23-24

# Best estimates of average time intensity and time-discrimination

  # beta estimates (Eq. 23)

  beta0 = colMeans(log(d[,5:29]),na.rm=TRUE)

  # alpha estimates (Eq. 24)

  alpha0 = sqrt(1/colMeans((log(d[,5:29]) - matrix(beta0,93,25,byrow=T))^2,na.rm=TRUE))

  
  # hyper parameters for inverse gamma 
  # Levy, Bayesian Psychometric Analysis, page 83, inverse gamma prior distribution
  
  N= 93
  a = mean(alpha0)

  v0 = 93/2
  v0
  v0*a

    # alpha ~ InvGamm(46.5, 52.8)
 
   require(invgamma)

   test.gamma = rinvgamma(10000,46.5,52.8)
    hist(test.gamma)
    mean(test.gamma)
    
    test.gamma = rgamma(10000,46.5,52.8)
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
    
    

################################################################################
# STAN model syntax for the DG-LNRT model
    
'    
    data{
      int <lower=1> I;                           // number of individuals              
      int <lower=1> J;                           // number of items
      int <lower=1> n_obs;                       // total number of observations (I x J) excluding missing responses
      int <lower=1> ind_person_obs[n_obs];       // person position indicator for an observed response
      int <lower=1> ind_item_obs[n_obs];         // item position indicator for an observed response
      int <lower=0,upper=1> i_status_obs[n_obs]; // item status indicator (1: disclosed, 0: nondisclosed)
      real Y[n_obs];                             // vector of response times
    }
    
    parameters {
      vector[J] beta;               // vector of time-intensity parameters             
      vector <lower=0> [J]  alpha;  // vector of time-discrimination parameters 
      vector[I] tau_t;              // vector of latent speed parameters for non-disclosed items
      vector[I] tau_c;              // vector of latent speed parameters for disclosed items
      real mu1;                     // hyperparameter for the mean of betas
      real<lower=0> sigma1;         // hyperparameter for the sd of betas
      real<lower=0> sigma_t;        // standard deviation of tau_t
      real<lower=0> sigma_c;        // standard deviation of tau_c
    }
    
    transformed parameters {
      vector[I] T;                  // vector of T for posterior estimate of item preknowledge
                                    // for each sampling, compare tau_t and tau_c
                                    // if tau_c > tau_t, then T=1, 0 otherwise
      
      for (i in 1:n_obs) {
        if(tau_t[ind_person_obs[i]]>tau_c[ind_person_obs[i]])
          
          T[ind_person_obs[i]] = 0;
        
        else 
          
          T[ind_person_obs[i]] = 1;
      }
    }
    
    
    model{
      
      sigma_t ~ exponential(1);      // Prior for sigma_t
      sigma_c ~ exponential(1);      // Prior for sigma_c
      
      tau_t    ~ normal(0,sigma_t);  // Prior for tau_t, mean is fixed to 0
      tau_c    ~ normal(0,sigma_c);  // Prior for tau_c, mean is fixed to 0
      
      mu1      ~ normal(3.76,1);     // hyperPrior for mean of betas
                                     // 3.76 is estimated from data see above under Prior
                                     
      sigma1   ~ exponential(1);     // hyperPrior for sigma1
      
      beta     ~ normal(mu1,sigma1); // prior for beta with hyperparameters mu1 and sigma1
      
      alpha    ~ inv_gamma(46.5,52.8); // prior for alpha
                                       // 46.5 and 52.8 are estimated from data see above under Prior
                                     
      
      // sample a response from a normal distribution, N(beta - tau, 1/alpha)
      // for disclosed items use tau_c
      // for nondisclosed items use tau_t
      
      for (i in 1:n_obs) {
        
        real p_t = beta[ind_item_obs[i]]-tau_t[ind_person_obs[i]];
        real p_c = beta[ind_item_obs[i]]-tau_c[ind_person_obs[i]];
        
        real p = (p_t^(1-T[ind_person_obs[i]]))*
          (((1-i_status_obs[i])*p_t + 
              (i_status_obs[i])*p_c)^T[ind_person_obs[i]]);
        
        Y[i] ~ normal(p,1/(alpha[ind_item_obs[i]]));
      }
    }
'
    
###############################################################################

# Read the Stan Model Syntax 
    
mod <- cmdstan_model(here('dglnrt2.stan'))

# Fit the model using cmdstan
    
fit <- mod$sample(data = data_rt,
                  seed = 1234,
                  chains = 4,
                  parallel_chains = 4,
                  iter_warmup   = 1000,
                  iter_sampling = 4000,
                  refresh = 100
                  )

# Save the output object

fit$cmdstan_summary()

stanfit <- rstan::read_stan_csv(fit$output_files())

save.image('dglnrt_v1/results with original priors.RData')

################################################################################
#                               OUTPUT ANALYSIS

# Load the pre-saved model output

load(here('data/dglnrt_v1/results with original priors.RData'))

# Summary of all parameters

params <- summary(stanfit, 
                  pars = c("mu1","sigma1","sigma_t","sigma_c",
                           'beta','alpha','tau_t','tau_c','T'), 
                  probs = c(0.025, 0.975))$summary

# Extract the betas

betas <- summary(stanfit, pars = c("beta"), probs = c(0.025, 0.975))$summary
betas 
describe(betas[,1])

# Extract the alphas

alphas <- summary(stanfit, pars = c("alpha"), probs = c(0.025, 0.975))$summary
alphas
describe(alphas[,1])

# Extract the tau_t

tau_t <- summary(stanfit, pars = c("tau_t"), probs = c(0.025, 0.975))$summary
tau_t
describe(tau_t[,1])

# Extract the tau_c

tau_c <- summary(stanfit, pars = c("tau_c"), probs = c(0.025, 0.975))$summary
tau_c
describe(tau_c[,1])

# Extract Ts

Ts <- summary(stanfit, pars = c("T"), probs = c(0.025, 0.975))$summary
Ts

# For a given threshold, check the number of identified individuals in the data 
# and compare it to their known status of item preknowledge

th = .99

table(ifelse(Ts[,1]>th,1,0),d.sub$COND)

hist(Ts[,1])


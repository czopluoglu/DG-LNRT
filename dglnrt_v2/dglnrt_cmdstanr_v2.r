require(cmdstanr)
require(here)
require(psych)
require(rstan)
require(pROC)
################################################################################

# This fits the model to the real dataset by Cizek and Wollack

################################################################################

# Import the data in the long format

d.long <- read.csv(here('data/rt_long.csv'))

# Remove missing data in the long format.
# If someone didn't respond to an item, only that particular item was removed
# the rest of the data for other items for the person remained.

d.long <- na.omit(d.long)

# Add a variable for a numeric item number based on item label

d.long$item <- as.numeric(substring(d.long$Item,6))

# a vector for the list of uncommon items in Form 2

icode = c(201,204,205,206,208,210,211,214,216,218,220,222,223,224,226,227,228,
          229,230,231,232,234,236,237,239,241,242,243,244,247,248,249,253,255,
          257,258,261,263,264,265,267,268,269,272,279,280,281,282,283,286,287,
          288,292,296,297,298,302,303,305,307,308,314,316,318,322,323,327,329,
          331,338,340,341,343,345,346,347,348,350,353,359,366,367,368)

icode2 = 171:253


# Change the item number for uncommon items in a way that the item numbers go
# from 1 to 253

for(i in 1:83){
  d.long[which(d.long$item == icode[i]),]$item = icode2[i]
}

unique(d.long$item)

  # Quality check to make sure things were properly changed

    tab <- as.data.frame(table(d.long$Item,d.long$item))
    tab = tab[which(tab[,3]!=0),]
    tab

# Add the id variable based on unique EID
    
d.long$id <- NA

ind = which(substring(d.long$EID,1,3)=='e10')
d.long[ind,]$id = as.numeric(substring(d.long[ind,]$EID,4,7))

ind = which(substring(d.long$EID,1,3)=='e20')
d.long[ind,]$id = as.numeric(substring(d.long[ind,]$EID,4,7))+1636

  # Quality check to make sure things were properly changed

    tab <- as.data.frame(table(d.long$EID,d.long$id))
    tab = tab[which(tab[,3]!=0),]
    tab

    unique(d.long$id)
########################################################################

# Add a variable for item compromised status
    
d.long$i.status <- ifelse(d.long$i_flag=='Flagged',1,0)

# Data object for Stan run
    
data_rt <- list(J              = 253,
                I              = 3280,
                n_obs          = length(d.long$RT),
                ind_person_obs = d.long$id,
                ind_item_obs   = d.long$item,
                i_status_obs   = d.long$i.status,
                Y              = d.long$RT)

################################################################################

# Stan Model Syntax

'
data{
  int <lower=1> I;                           // number of individuals                               
  int <lower=1> J;                           // number of items
  int <lower=1> n_obs;                       // number of observations (I x J)
  int <lower=1> ind_person_obs[n_obs];       // person location indicator
  int <lower=1> ind_item_obs[n_obs];         // item location indicator
  int <lower=0,upper=1> i_status_obs[n_obs]; // item compromised status indicator
  real Y[n_obs];                             // observed response vector
}

parameters {
  vector[J] beta;                           // time-intensity parameters
  vector <lower=0> [J]  alpha;              // time-discrimination parameters
  vector[I] tau_t;                          // latent speeed parameters for uncompromised items
  vector[I] tau_c;                          // latent speed parameters for compromised items 
  real mu1;                                 // hyperparameter for mean of beta parameters
  real<lower=0> sigma1;                     // hyperparameter for sd of beta parameters
  real<lower=0> sigma_t;                    // standard deviation of tau_t
  real<lower=0> sigma_c;                    // standard deviation of tau_c
}

transformed parameters {
  vector[I] T;               // vector of T for posterior estimate of item preknowledge
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
  
  sigma_t ~ exponential(1);         // Prior for sigma_t
  sigma_c ~ exponential(1);         // Prior for sigma_c
  
  tau_t    ~ normal(0,sigma_t);     // Prior for tau_t, mean is fixed to 0
  tau_c    ~ normal(0,sigma_c);     // Prior for tau_c, mean is fixed to 0
  
  mu1      ~ normal(3.98,1);        // hyperPrior for mean of betas
                                    // 3.98 is estimated from data see above under Prior
                                    
  sigma1   ~ exponential(1);        // hyperPrior for sigma1
  beta     ~ normal(mu1,sigma1);    // prior for beta with hyperparameters mu1 and sigma1
  
  alpha    ~ inv_gamma(800,1550);   // prior for alpha
                                    // 800 and 1550 are estimated from data see above under Prior
                                     
  
  // sample a response from a normal distribution, N(beta - tau, 1/alpha)
  // for disclosed items use tau_c
  // for nondisclosed items use tau_t
      
  for (i in 1:n_obs) {
    
    real p_t = beta[ind_item_obs[i]]-tau_t[ind_person_obs[i]];
    real p_c = beta[ind_item_obs[i]]-tau_c[ind_person_obs[i]];
    
    real p = (p_t^(1-T[ind_person_obs[i]]))*
      (((1-i_status_obs[i])*p_t + 
          (i_status_obs[i])*p_c)^T[ind_person_obs[i]]);
    
    Y[i] ~ normal(p,1/(alpha[ind_item_obs[i]]^2));
  }
}

'

################################################################################

# PRIORS


# van der Linder 2006
# page 191, Eq 23-24

# Best estimates of average time intensity and time-discrimination

  # beta estimates (Eq. 23)

  beta0 = describeBy(d.long$RT,d.long$item,mat=TRUE)$mean
  
  # alpha estimates (Eq. 24)
  
  alpha0 = c()
  for(i in 1:253){
    alpha0[i] = sqrt(1/mean((d.long[which(d.long$item==i),]$RT - beta0[i])^2))
  }

  mean(beta0)
  mean(alpha0)
  
  # hyper parameters for inverse gamma 
  # Levy, Bayesian Psychometric Analysis, page 83, inverse gamma prior distribution
  
  N= 1600
  a = mean(alpha0)

  v0 = N/2
  v0
  v0*a

  # alpha ~ InvGamm(800,1550)

require(invgamma)

test.gamma = rinvgamma(10000,v0,v0*a)
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

# Read the Stan Model Syntax 

mod <- cmdstan_model(here('dglnrt_v2.stan'))

# Fit the model using cmdstan

fit <- mod$sample(
  data = data_rt,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_warmup   = 500,
  iter_sampling = 2000,
  refresh = 10,
  adapt_delta = 0.99,
  max_treedepth = 15
)


# Save the output object

fit$cmdstan_summary()

stanfit <- rstan::read_stan_csv(fit$output_files())

save.image(here('dglnrt_v2/results.RData'))


################################################################################
#                               OUTPUT ANALYSIS

# Load the pre-saved model output

load(here("data/dglnrt_v2/results.RData"))

# Extract Ts, posterior probability estimate of item preknowledge

T <- summary(stanfit, pars = c("T"), probs = c(0.025, 0.975))$summary
describe(T[,1])

# Create a vector for person status (flagged vs. not flagged) to compare with Ts

flagged <- c()
for(i in 1:3280){
  flagged[i] <- unique(d.long[which(d.long$id==i),]$p_flag)
}

table(flagged)


# For a given threshold, check the number of identified individuals in the data 
# and compare it to their known status of item preknowledge

comp <- cbind(flagged,T[,1])

th = 0.99

tab <- table(comp[,1],(T[,1]>th)*1)

round(tab[2,2]/sum(tab[2,]),3)  # tpr
round(tab[1,2]/sum(tab[1,]),3)  # fpr
round(tab[2,2]/sum(tab[,2]),3)  # pr


################################################################################
summary(stanfit, pars = c("mu1","sigma1","sigma_t","sigma_c"), 
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
























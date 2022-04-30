require(cmdstanr)
require(here)
require(psych)
require(rstan)
require(pROC)
################################################################################

# This fits the model to the real dataset (Form A) by Cizek and Wollack (2017)

################################################################################

# Import the data in the long format

d.long <- read.csv(here('data/rt_long.csv'))

# Remove the examinees from Form 2

d.long <- d.long[-grep('e20',d.long$EID),]

# Remove missing data in the long format.
# If someone didn't respond to an item, only that particular item was removed
# the rest of the data for other items for the person remained.

d.long <- na.omit(d.long)

# Add a variable for a numeric item number based on item label

d.long$item <- as.numeric(substring(d.long$Item,6))

# Check sample size per item

    tab <- as.data.frame(table(d.long$Item,d.long$item))
    tab = tab[which(tab[,3]!=0),]
    tab

# Check flagged items
    
    tab <- as.data.frame(table(d.long$Item,d.long$i_flag))
    tab = tab[which(tab[,3]!=0),]
    table(tab$Var2)
    
# Check flagged individuals
    
    tab <- as.data.frame(table(d.long$EID,d.long$p_flag))
    tab = tab[which(tab[,3]!=0),]
    table(tab$Var2)
    
# Add the id variable based on unique EID
    
d.long$id <- NA

ind = which(substring(d.long$EID,1,3)=='e10')
d.long[ind,]$id = as.numeric(substring(d.long[ind,]$EID,4,7))


  # Quality check to make sure things were properly changed

    tab <- as.data.frame(table(d.long$EID,d.long$id))
    tab = tab[which(tab[,3]!=0),]
    tab

    unique(d.long$id)
    
# There are 1636 examinees, 170 items
# 64 items flagged
# 46 individuals flagged
    
    
########################################################################
# Descriptive Statistics
    
# Flagged Items - Unflagged Examinees
    
flagged <- sort(unique(d.long[d.long$i_flag=='Flagged',]$item))
    
desc.rt <- matrix(nrow=length(flagged),ncol=4)

desc.rt[,1] <- flagged

for(i in 1:nrow(desc.rt)){
  
  temp <- exp(d.long[d.long$item==flagged[i] & d.long$p_flag==0,]$RT)
  
  desc.rt[i,2] <- length(temp)
  desc.rt[i,3] <- mean(temp)
  desc.rt[i,4] <- sd(temp)
}

rt10 <- desc.rt[,3]

# Flagged Items - Flagged Examinees


flagged <- sort(unique(d.long[d.long$i_flag=='Flagged',]$item))

desc.rt <- matrix(nrow=length(flagged),ncol=4)

desc.rt[,1] <- flagged

for(i in 1:nrow(desc.rt)){
  
  temp <- exp(d.long[d.long$item==flagged[i] & d.long$p_flag==1,]$RT)
  
  desc.rt[i,2] <- length(temp)
  desc.rt[i,3] <- mean(temp)
  desc.rt[i,4] <- sd(temp)
}

rt11 <- desc.rt[,3]


plot(rt10,rt11,
     ylim=c(10,120),
     xlim=c(10,120),
     xlab = 'Response Time in Seconds by Unflagged Examinees',
     ylab = 'Response Time in Seconds by Flagged Examinees')

abline(0,1,lty=2,col='gray')


# Unflagged Items - Unflagged Examinees

unflagged <- sort(unique(d.long[d.long$i_flag=='Unflagged',]$item))

desc.rt <- matrix(nrow=length(unflagged),ncol=4)

desc.rt[,1] <- unflagged

for(i in 1:nrow(desc.rt)){
  
  temp <- exp(d.long[d.long$item==unflagged[i] & d.long$p_flag==0,]$RT)
  
  desc.rt[i,2] <- length(temp)
  desc.rt[i,3] <- mean(temp)
  desc.rt[i,4] <- sd(temp)
}  


rt00 <- desc.rt[,3]



# Unflagged Items - Flagged Examinees

unflagged <- sort(unique(d.long[d.long$i_flag=='Unflagged',]$item))

desc.rt <- matrix(nrow=length(unflagged),ncol=4)

desc.rt[,1] <- unflagged

for(i in 1:nrow(desc.rt)){
  
  temp <- exp(d.long[d.long$item==unflagged[i] & d.long$p_flag==1,]$RT)
  
  desc.rt[i,2] <- length(temp)
  desc.rt[i,3] <- mean(temp)
  desc.rt[i,4] <- sd(temp)
}  


rt01 <- desc.rt[,3]


plot(rt00,rt01,
     ylim=c(10,120),
     xlim=c(10,120),
     xlab = 'Response Time in Seconds by Unflagged Examinees',
     ylab = 'Response Time in Seconds by Flagged Examinees',
     main = 'Unflagged Items')

abline(0,1,lty=2,col='gray')



mean(rt10)
mean(rt11)


mean(rt00)
mean(rt01)



########################################################################

# Add a variable for item compromised status
    
d.long$i.status <- ifelse(d.long$i_flag=='Flagged',1,0)

length(unique(d.long$id))

length(unique(d.long$item))




# Data object for Stan run
    
data_rt <- list(J              = 170,
                I              = 1636,
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
  for(i in 1:170){
    alpha0[i] = sqrt(1/mean((d.long[which(d.long$item==i),]$RT - beta0[i])^2))
  }

  mean(beta0)
  mean(alpha0)
  
  # hyper parameters for inverse gamma 
  # Levy, Bayesian Psychometric Analysis, page 83, inverse gamma prior distribution
  
  N= 1636
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

mod <- cmdstan_model(here('dglnrt_v2/dglnrt_v2.stan'))

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

save.image(here('data/dglnrt_v2/results.RData'))


################################################################################
#                               OUTPUT ANALYSIS

# Load the pre-saved model output

load(here("data/dglnrt_v2/results.RData"))

# Extract Ts, posterior probability estimate of item preknowledge

T <- summary(stanfit, pars = c("T"), probs = c(0.025, 0.975))$summary
describe(T[,1])


# Create a vector for person status (flagged vs. not flagged) to compare with Ts

flagged <- c()
for(i in 1:1636){
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
describe(alphas[,1]^2)

alpha.chain <- extract(stanfit)$alpha

alphas2 <- alpha.chain^2


# Extract the tau_t

tau_t <- summary(stanfit, pars = c("tau_t"), probs = c(0.025, 0.975))$summary
tau_t
describe(tau_t[,1])


# Extract the tau_c

tau_c <- summary(stanfit, pars = c("tau_c"), probs = c(0.025, 0.975))$summary
tau_c
describe(tau_c[,1])



hist(c(T[,7],betas[,7],alphas[,7],tau_t[,7],tau_c[,7]))

mean(c(T[,7],betas[,7],alphas[,7],tau_t[,7],tau_c[,7]),na.rm=TRUE)
min(c(T[,7],betas[,7],alphas[,7],tau_t[,7],tau_c[,7]),na.rm=TRUE)
max(c(T[,7],betas[,7],alphas[,7],tau_t[,7],tau_c[,7]),na.rm=TRUE)


describe(betas[,1])

describe(alphas[,1]^2)
































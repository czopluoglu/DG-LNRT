
require(cmdstanr)
require(here)
require(psych)
########################

# This fits the model to the real dataset by Toton and Maynes (2019)

#########################

d <- read.csv(here('data/uva_rt.csv'),na.strings = '#N/A')


d.sub <- d[,c("ï..ID","COND","Q1RT","Q2RT","Q3RT","Q4RT","Q5RT","Q6RT","Q7RT","Q8RT","Q9RT","Q10RT",
              "Q11RT","Q12RT","Q13RT","Q14RT","Q15RT","Q16RT","Q17RT","Q18RT","Q19RT","Q20RT",
              "Q21RT","Q22RT","Q23RT","Q24RT","Q25RT")]

colnames(d.sub)[1] = 'ID'

d.sub$ID <- 1:93

d.long <- reshape(data        = d.sub,
                  idvar       = "ID",
                  varying     = list(colnames(d.sub)[3:27]),
                  timevar     = "Item",
                  times       = 1:25,
                  v.names      = "RT",
                  direction   = "long")

d.long$i.status <- ifelse(d.long$Item%%2==0,1,0)

d.long <- na.omit(d.long)

d.long$logRT <- log(d.long$RT)

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

  beta0 = colMeans(log(d[,5:29]),na.rm=TRUE)

  alpha0 = sqrt(1/colMeans((log(d[,5:29]) - matrix(beta0,93,25,byrow=T))^2,na.rm=TRUE))

  mean(beta0)
  mean(alpha0)

  
  N= 93
  a = mean(alpha0)

  # hyper parameters for inverse gamma 
  # Levy, Bayesian Psychometric Analysis, page 83, inverse gamma prior distribution
  
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
    
    

##########################################################

mod <- cmdstan_model(here('dglnrt2.stan'))

fit <- mod$sample(
  data = data_rt,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_warmup   = 1000,
  iter_sampling = 4000,
  refresh = 100
)

fit$cmdstan_summary()

stanfit <- rstan::read_stan_csv(fit$output_files())

save.image('dglnrt_v1/results with original priors.RData')

####################################################

summary(stanfit, pars = c("mu1","sigma1","sigma_t","sigma_c"), 
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

tau_c <- summary(stanfit, pars = c("tau_c"), probs = c(0.025, 0.975))$summary
tau_c
describe(tau_c[,1])

Ts <- summary(stanfit, pars = c("T"), probs = c(0.025, 0.975))$summary
Ts

th = .95

table(ifelse(Ts[,1]>th,1,0),d.sub$COND)

hist(Ts[,1])



####################################################















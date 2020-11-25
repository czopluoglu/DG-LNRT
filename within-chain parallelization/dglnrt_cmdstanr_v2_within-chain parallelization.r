require(cmdstanr)
require(here)
require(psych)
########################

# This fits the model to the real dataset by Cizek and Wollack

#########################


d.long <- read.csv(here('data/rt_long.csv'))

d.long <- na.omit(d.long)

d.long$item <- as.numeric(substring(d.long$Item,6))

icode = c(201,204,205,206,208,210,211,214,216,218,220,222,223,224,226,227,228,
          229,230,231,232,234,236,237,239,241,242,243,244,247,248,249,253,255,
          257,258,261,263,264,265,267,268,269,272,279,280,281,282,283,286,287,
          288,292,296,297,298,302,303,305,307,308,314,316,318,322,323,327,329,
          331,338,340,341,343,345,346,347,348,350,353,359,366,367,368)

icode2 = 171:253


for(i in 1:83){
  d.long[which(d.long$item == icode[i]),]$item = icode2[i]
}

unique(d.long$item)

tab <- as.data.frame(table(d.long$Item,d.long$item))
tab = tab[which(tab[,3]!=0),]



d.long$id <- NA

ind = which(substring(d.long$EID,1,3)=='e10')
d.long[ind,]$id = as.numeric(substring(d.long[ind,]$EID,4,7))

ind = which(substring(d.long$EID,1,3)=='e20')
d.long[ind,]$id = as.numeric(substring(d.long[ind,]$EID,4,7))+1636

tab <- as.data.frame(table(d.long$EID,d.long$id))
tab = tab[which(tab[,3]!=0),]


unique(d.long$id)
########################################################################

d.long$i.status <- ifelse(d.long$i_flag=='Flagged',1,0)

data_rt <- list(J              = 253,
                I              = 3280,
                n_obs          = length(d.long$RT),
                ind_person_obs = d.long$id,
                ind_item_obs   = d.long$item,
                i_status_obs   = d.long$i.status,
                i_status_obs2  = 1- d.long$i.status,
                Y              = d.long$RT)

########################################################################
##########################################################

# PRIORS


# van der Linder 2006
# page 191, Eq 23-24

# Best estimates of average time intensity and time-discrimination

beta0 = describeBy(d.long$RT,d.long$item,mat=TRUE)$mean
  
alpha0 = c()
for(i in 1:253){
  alpha0[i] = sqrt(1/mean((d.long[which(d.long$item==i),]$RT - beta0[i])^2))
}

mean(beta0)
mean(alpha0)


N= 1600
a = mean(alpha0)

# hyper parameters for inverse gamma 
# Levy, Bayesian Psychometric Analysis, page 83, inverse gamma prior distribution

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



##########################################################

mod <- cmdstan_model(here('within-chain parallelization/dglnrt_cmdstanr_v2_within-chain parallelization.stan'))

fit <- mod$sample(
  data = data_rt,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 5,
  iter_warmup   = 100,
  iter_sampling = 250,
  refresh = 10,
  adapt_delta = 0.99
)

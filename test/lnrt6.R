
require(rstan)
require(cmdstanr)
require(here)
######################################################

b = rnorm(10,4,0.5)
a = runif(10,1,2)
t <- rnorm(150,0,.5)

d <- matrix(nrow=150,ncol=10)

for(i in 1:150){
  for(j in 1:10){
    d[i,j]=rnorm(1,b[j]-t[i],1/a[j])
  }
}

d <- as.data.frame(d)

d$id <- 1:150

d.long <- reshape(d,
                  varying=list(colnames(d)[1:10]),
                  idvar = 'id',
                  direction = 'long',
                  v.names = 'resp')

#######################################################

data_rt <- list(I              = 150,
                J              = 10,
                n_obs          = 1500,
                p              = d.long$id,
                it             = d.long$time,
                Y              = d.long$resp,
                grainsize      = 1)


mod <- cmdstan_model(here('test/lnrt6.stan'),
                     cpp_options = list(stan_threads = TRUE))

fit <- mod$sample(
  data = data_rt,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 2,
  iter_warmup   = 500,
  iter_sampling = 2000,
  refresh = 10,
  adapt_delta = 0.99,
  max_treedepth = 15
)

fit$cmdstan_summary()

stanfit <- rstan::read_stan_csv(fit$output_files())


summary(stanfit, pars = c("mua","sigmaa","mub","sigmab","sigmat"), 
        probs = c(0.025, 0.975))$summary

betas <- summary(stanfit, pars = c("b"), probs = c(0.025, 0.975))$summary
betas 
cbind(b,betas[,1])

alphas <- summary(stanfit, pars = c("a"), probs = c(0.025, 0.975))$summary
alphas
cbind(a,alphas[,1])

tau <- summary(stanfit, pars = c("t"), probs = c(0.025, 0.975))$summary
plot(cbind(t,tau[,1]))
cor(cbind(t,tau[,1]))


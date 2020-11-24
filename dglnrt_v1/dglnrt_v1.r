
require(here)
require(rstan)

########################

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

code_rt <- '

  data{
    int <lower=1> I;                      //  number of individuals
    int <lower=1> J;                      //  number of items
    int <lower=1> n_obs;                  //  number of observed responses
    int <lower=1> ind_person_obs[n_obs];  //  person indicator of a response
    int <lower=1> ind_item_obs[n_obs];    //  item indicator of a response
    int <lower=0,upper=1> i_status_obs[n_obs]; //  item status indicator of a response 
    real Y[n_obs];                       //  vector of response times
  }

  parameters {
    vector[J] beta;             //item specific time intensity parameter
    vector <lower=0> [J]  alpha;  // item specific residual standard deviation,
    vector[I] tau_t;            // person specific latent working speed
    vector[I] tau_c;            // person specific latent working speed
    real mu1;
    real<lower=0> sigma1;
  }

  transformed parameters {
    vector[I] T;
  
    for (i in 1:n_obs) {
      if(tau_t[ind_person_obs[i]]>tau_c[ind_person_obs[i]])
        
        T[ind_person_obs[i]] = 0;
      
      else 
        
        T[ind_person_obs[i]] = 1;
    }
  }


  model{
  
    tau_t    ~ normal(0,1);
    tau_c    ~ normal(0,1);
    
    mu1      ~ normal(0,1);
    sigma1   ~ exponential(1);
    
    beta     ~ normal(mu1,sigma1);
    alpha    ~ inv_gamma(5,5);         
    
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

####################################################



rt_stan <- stan(model_code = code_rt, 
                data = data_rt, 
                iter =5000,
                chains = 3,
                control = list(max_treedepth = 15),
                cores = getOption("mc.cores", 4L))



sum(get_elapsed_time(rt_stan))/(3600)

####################################################


betas <- summary(rt_stan, pars = c("beta"), probs = c(0.025, 0.975))$summary
betas 

alphas <- summary(rt_stan, pars = c("alpha"), probs = c(0.025, 0.975))$summary
alphas

tau_t <- summary(rt_stan, pars = c("tau_t"), probs = c(0.025, 0.975))$summary
tau_t

tau_c <- summary(rt_stan, pars = c("tau_c"), probs = c(0.025, 0.975))$summary
tau_c

Ts <- summary(rt_stan, pars = c("T"), probs = c(0.025, 0.975))$summary
Ts

th = .5

table(ifelse(Ts[,1]>th,1,0),d.sub$COND)




















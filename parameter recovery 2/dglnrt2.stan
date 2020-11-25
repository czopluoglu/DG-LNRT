data{
  int <lower=1> I;                      
  int <lower=1> J;                      
  int <lower=1> n_obs;                  
  int <lower=1> ind_person_obs[n_obs];  
  int <lower=1> ind_item_obs[n_obs];    
  real Y[n_obs];                       
}

parameters {
  vector[J] beta;             
  vector <lower=0> [J]  alpha;  
  vector[I] tau_t;            
  real mu1;
  real<lower=0> sigma1;
  real<lower=0> sigma_t;
}

model{
  
  sigma_t ~ exponential(1);
    tau_t    ~ normal(0,sigma_t);
  
  mu1      ~ normal(4.20,1);
  sigma1   ~ exponential(1);
    beta     ~ normal(mu1,sigma1);

    alpha    ~ inv_gamma(46.5,71.96);
  
  for (i in 1:n_obs) {
    Y[i] ~ normal(beta[ind_item_obs[i]]-tau_t[ind_person_obs[i]],1/(alpha[ind_item_obs[i]]));
  }
}
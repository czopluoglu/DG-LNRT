data{
  int <lower=1> I;                      
  int <lower=1> J;                      
  int <lower=1> n_obs;                  
  int <lower=1> ind_person_obs[n_obs];  
  int <lower=1> ind_item_obs[n_obs];    
  int <lower=0,upper=1> i_status_obs[n_obs];
  real Y[n_obs];                       
}

parameters {
  vector[J] beta;             
  vector <lower=0> [J]  alpha;  
  vector[I] tau_t;            
  vector[I] tau_c;            
  real mu1;
  real<lower=0> sigma1;
  real<lower=0> sigma_t;
  real<lower=0> sigma_c;
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
  
  sigma_t ~ exponential(1);
  sigma_c ~ exponential(1);
  
  tau_t    ~ normal(0,sigma_t);
  tau_c    ~ normal(0,sigma_c);
  
  mu1      ~ normal(3.97,1);
  sigma1   ~ exponential(1);
    beta     ~ normal(mu1,sigma1);

    alpha    ~ inv_gamma(818,1571);
  
  for (i in 1:n_obs) {
    
    real p_t = beta[ind_item_obs[i]]-tau_t[ind_person_obs[i]];
    real p_c = beta[ind_item_obs[i]]-tau_c[ind_person_obs[i]];
    
    real p = (p_t^(1-T[ind_person_obs[i]]))*
      (((1-i_status_obs[i])*p_t + 
          (i_status_obs[i])*p_c)^T[ind_person_obs[i]]);
    
    Y[i] ~ normal(p,1/(alpha[ind_item_obs[i]]));
  }
}
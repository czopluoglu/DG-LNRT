functions {
  real partial_sum(real[] slice_Y,
                   int start, int end,
                   int[] ind_person_obs,
                   int[] ind_item_obs,
                   int[] i_status_obs,
                   vector beta,
                   vector alpha,
                   vector tau_t,
                   vector tau_c,
                   vector T) {
    
    return normal_lpdf(slice_Y | ((beta[ind_item_obs[start:end]]-tau_t[ind_person_obs[start:end]]).^(1 .- T[ind_person_obs[start:end]]))*(((1 .- i_status_obs[start:end])*(beta[ind_item_obs[start:end]]-tau_t[ind_person_obs[start:end]]) + (i_status_obs[start:end])*beta[ind_item_obs[start:end]]-tau_c[ind_person_obs[start:end]])^T[ind_person_obs[start:end]]), 1);
    }
  }


data{
  int <lower=1> I;                      
  int <lower=1> J;                      
  int <lower=1> n_obs;                  
  int <lower=1> ind_person_obs[n_obs];  
  int <lower=1> ind_item_obs[n_obs];    
  int <lower=0,upper=1> i_status_obs[n_obs];
  real Y[n_obs];      
  int<lower=1> grainsize;
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
  
  mu1      ~ normal(3.98,1);
  sigma1   ~ exponential(1);
  beta     ~ normal(mu1,sigma1);
  
  alpha    ~ inv_gamma(800,1550);
  
  target += reduce_sum(partial_sum, Y , grainsize, ind_person_obs, ind_item_obs,
                       i_status_obs, beta, alpha, tau_t, tau_c, T);
}








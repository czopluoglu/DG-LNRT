data{
  int <lower=1> I;                           // number of examinees
  int <lower=1> J;                           // number of items
  int <lower=1> n_obs;                       // number of observations
  int <lower=1> ind_person_obs[n_obs];       // examinee position index
  int <lower=1> ind_item_obs[n_obs];         // item position index
  int <lower=0,upper=1> i_status_obs[n_obs]; // item status index
  int <lower=0,upper=1> Y[n_obs];            // observed binary responses
}

parameters {
 real mu_thetat;             // mean for theta_t
 real<lower=0> sigma_thetat; // sd for theta_t
  
 real mu_thetac;             // mean for theta_c
 real<lower=0> sigma_thetac; // sd for theta_c

 vector[I] theta_t;           // vector of true theta parameters
 vector[I] theta_c;           // vector of cheating theta parameters
 vector[J] b;                // vector of item difficulty parameters

}


transformed parameters {
  vector[I] T;
  
  for (i in 1:n_obs) {
    if(theta_t[ind_person_obs[i]]>theta_c[ind_person_obs[i]])
      
      T[ind_person_obs[i]] = 0;
    
    else 
      
      T[ind_person_obs[i]] = 1;
  }
}


model{
  
  sigma_thetat ~ exponential(1);
  sigma_thetac ~ exponential(1);
  mu_thetat    ~ normal(0,10);
  mu_thetac    ~ normal(0,10);

  theta_t    ~ normal(mu_thetat,sigma_thetat);
  theta_c    ~ normal(mu_thetac,sigma_thetac);

  b          ~ normal(0,1);

  for (i in 1:n_obs) {
    
    real p_t = inv_logit(theta_t[ind_person_obs[i]] - b[ind_item_obs[i]]);
    real p_c = inv_logit(theta_c[ind_person_obs[i]] - b[ind_item_obs[i]]);
    
    real p = (p_t^(1-T[ind_person_obs[i]]))*
      (((1-i_status_obs[i])*p_t + 
          (i_status_obs[i])*p_c)^T[ind_person_obs[i]]);
    
    Y[i] ~ bernoulli(p);
  }
}


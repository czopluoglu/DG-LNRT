data{
  int <lower=1> I;                      
  int <lower=1> J;                      
  int <lower=1> n_obs;                  
  int <lower=1> p[n_obs];  
  int <lower=1> it[n_obs];  
  real Y[n_obs];                       
}

parameters {
  vector[J] b;             
  vector <lower=0> [J]  a;  
  vector[I] t;            
  real<lower=0> sigmat;
  real mua;
  real<lower=0> sigmaa;
  real mub;
  real<lower=0> sigmab;
}


model{
  
  mub      ~ normal(0,5);
  sigmab   ~ exponential(1);
  b     ~ normal(mub,sigmab);
  
  mua      ~ normal(0,5);
  sigmaa   ~ exponential(1);
  a     ~ normal(mua,sigmaa);
  
  sigmat   ~ exponential(1);
  t     ~ normal(0,sigmat);
  
  target += normal_lpdf(Y | b[it] - t[p], 1 ./ a[it]);
  
}
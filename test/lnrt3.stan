data{
  int <lower=1> I;                      
  int <lower=1> J;            
  matrix [I,J] Y;                       
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
  
  Y[,1] ~ normal(b[1] - t,1/a[1]);
  Y[,2] ~ normal(b[2] - t,1/a[2]);
  Y[,3] ~ normal(b[3] - t,1/a[3]);
  Y[,4] ~ normal(b[4] - t,1/a[4]);
  Y[,5] ~ normal(b[5] - t,1/a[5]);
  Y[,6] ~ normal(b[6] - t,1/a[6]);
  Y[,7] ~ normal(b[7] - t,1/a[7]);
  Y[,8] ~ normal(b[8] - t,1/a[8]);
  Y[,9] ~ normal(b[9] - t,1/a[9]);
  Y[,10] ~ normal(b[10] - t,1/a[10]);
  
}
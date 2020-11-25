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
  
  for (i in 1:I) {
    for(j in 1:J){
     Y[i,j] ~ normal(b[j] - t[i],1/a[j]); 
    }
     
  }

}
data{
  int <lower=1> I;                      
  int <lower=1> J;                      
  int <lower=1> n1;
  int <lower=1> n2;
  int <lower=1> n3;
  int <lower=1> n4;
  int <lower=1> n5;
  int <lower=1> n6;
  int <lower=1> n7;
  int <lower=1> n8;
  int <lower=1> n9;
  int <lower=1> n10;
  int <lower=1> n11;
  int <lower=1> n12;
  int <lower=1> n13;
  int <lower=1> n14;
  int <lower=1> n15;
  int <lower=1> p1[n1];
  int <lower=1> p2[n2];
  int <lower=1> p3[n3];
  int <lower=1> p4[n4];
  int <lower=1> p5[n5];
  int <lower=1> p6[n6];
  int <lower=1> p7[n7];
  int <lower=1> p8[n8];
  int <lower=1> p9[n9];
  int <lower=1> p10[n10];
  int <lower=1> p11[n11];
  int <lower=1> p12[n12];
  int <lower=1> p13[n13];
  int <lower=1> p14[n14];
  int <lower=1> p15[n15];
  real y1[n1];        
  real y2[n2];        
  real y3[n3];        
  real y4[n4];        
  real y5[n5];        
  real y6[n6];        
  real y7[n7];        
  real y8[n8];        
  real y9[n9];        
  real y10[n10];        
  real y11[n11];        
  real y12[n12];        
  real y13[n13];        
  real y14[n14];        
  real y15[n15];        
}

parameters {
  vector[I] beta;             
  vector <lower=0> [I]  alpha;  
  vector[J] tau_t;            
  vector[J] tau_c;            
  real mu1;
  real<lower=0> sigma1;
  real<lower=0> sigma_t;
  real<lower=0> sigma_c;
}

transformed parameters {
  vector[J] T;
  
  for (j in 1:J) {
    if(tau_t[j]>tau_c[j])
      
      T[j] = 0;
    
    else 
      
      T[j] = 1;
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
  
  vector[n1] m1 = ((beta[1] - tau_t[p1]) .^ (1 - T[p1])) .* ((beta[1] - tau_c[p1]) .^ T[p1]);
  y1 ~ normal(m1,1/alpha[1]);
 
  y2 ~ normal(beta[2]-tau_t[p2],1/alpha[2]);
  
  vector[n3] m3 = ((beta[3] - tau_t[p3]) .^ (1 - T[p3])) .* ((beta[3] - tau_c[p3]) .^ T[p3]);
  y3 ~ normal(m3,1/alpha[3]);
  
  vector[n4] m4 = ((beta[4] - tau_t[p4]) .^ (1 - T[p4])) .* ((beta[4] - tau_c[p4]) .^ T[p4]);
  y4 ~ normal(m4,1/alpha[4]);
  
  y5 ~ normal(beta[5]-tau_t[p5],1/alpha[5]);

  y6 ~ normal(beta[6]-tau_t[p6],1/alpha[6]);
  
  vector[n7] m7 = ((beta[7] - tau_t[p7]) .^ (1 - T[p7])) .* ((beta[7] - tau_c[p7]) .^ T[p7]);
  y7 ~ normal(m7,1/alpha[7]);
  
  vector[n8] m8 = ((beta[8] - tau_t[p8]) .^ (1 - T[p8])) .* ((beta[8] - tau_c[p8]) .^ T[p8]);
  y8 ~ normal(m8,1/alpha[8]);
  
  vector[n9] m9 = ((beta[9] - tau_t[p9]) .^ (1 - T[p9])) .* ((beta[9] - tau_c[p9]) .^ T[p9]);
  y9 ~ normal(m9,1/alpha[9]);
  
  vector[n10] m10 = ((beta[10] - tau_t[p10]) .^ (1 - T[p10])) .* ((beta[10] - tau_c[p10]) .^ T[p10]);
  y10 ~ normal(m10,1/alpha[10]);
  
  vector[n11] m11 = ((beta[11] - tau_t[p11]) .^ (1 - T[p11])) .* ((beta[11] - tau_c[p11]) .^ T[p11]);
  y11 ~ normal(m11,1/alpha[11]);
  
  vector[n12] m12 = ((beta[12] - tau_t[p12]) .^ (1 - T[p12])) .* ((beta[12] - tau_c[p12]) .^ T[p12]);
  y12 ~ normal(m12,1/alpha[12]);
  
  vector[n13] m13 = ((beta[13] - tau_t[p13]) .^ (1 - T[p13])) .* ((beta[13] - tau_c[p13]) .^ T[p13]);
  y13 ~ normal(m13,1/alpha[13]);
  
  y14 ~ normal(beta[14]-tau_t[p14],1/alpha[14]);
  
  y15 ~ normal(beta[15]-tau_t[p15],1/alpha[15]);
  
}


















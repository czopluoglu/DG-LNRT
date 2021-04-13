####################################################################
####################################################################
####################################################################

# A generic function to simulate DG-LNRT data with no item preknowledge

sim_dglnrt <- function() {
  
  # MODEL PARAMETERS
  
  beta  <- rnorm(25,4.15,0.39)
  alpha <- rnorm(25,1.36,0.30)
  
  # Tau for control group
  
  cor0 <- matrix(c(1,.84,.84,1),2,2)
  tau0 <- mvrnorm(33,c(0,0),cor2cov(cor0,c(0.18,0.17)))
  
  # Tau for experimental group (Items Disclosed)
  
  cor1 <- matrix(c(1,.56,.56,1),2,2)
  tau1 <- mvrnorm(30,c(0.05,1.34),cor2cov(cor1,c(0.55,0.47)))
  
  # Tau for experimental group (Items and Answers Disclosed)
  
  cor2 <- matrix(c(1,.43,.43,1),2,2)
  tau2 <- mvrnorm(30,c(-0.18,1.47),cor2cov(cor2,c(0.43,0.46)))
  
  # A vector for item status (0: not disclosed, 1:disclosed)
  
  C    <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  
  # RESPONSE TIME GENERATION
  
  # Note that the response time data generated is
  # already on the log scale
  
  # Control Group
  
  rt0 <- matrix(nrow = 33, ncol = 25)
  
  for (i in 1:33) {
    for (j in 1:25) {
      p_t = beta[j] - tau0[i, 1]
      p_c = beta[j] - tau0[i, 2]
      p   = p_t * (1 - C[j]) + p_c * C[j]
      rt0[i, j] = rnorm(1, p, 1 / alpha[j])
    }
  }
  
  # Experimental Group 1 (Item Disclosed)
  
  rt1 <- matrix(nrow = 30, ncol = 25)
  
  for (i in 1:30) {
    for (j in 1:25) {
      p_t = beta[j] - tau1[i, 1]
      p_c = beta[j] - tau1[i, 2]
      p   = p_t * (1 - C[j]) + p_c * C[j]
      rt1[i, j] = rnorm(1, p, 1 / alpha[j])
    }
  }
  
  # Experimental Group 2 (Item and Answers Disclosed)
  
  rt2 <- matrix(nrow = 30, ncol = 25)
  
  for (i in 1:30) {
    for (j in 1:25) {
      p_t = beta[j] - tau2[i, 1]
      p_c = beta[j] - tau2[i, 2]
      p   = p_t * (1 - C[j]) + p_c * C[j]
      rt2[i, j] = rnorm(1, p, 1 / alpha[j])
    }
  }
  
  # Combine the groups
  
  rt <- rbind(cbind(data.frame(exp(rt0)), gr = 1),
              cbind(data.frame(exp(rt1)), gr = 2),
              cbind(data.frame(exp(rt2)), gr = 3))
  
  
  return(list(
    rt = rt,
    b = beta,
    a = alpha,
    tau_t = c(tau0[, 1], tau1[, 1], tau2[,1]),
    tau_c = c(tau0[, 2], tau1[, 2], tau2[,2])
  ))
}



PPest <- function(alpha, beta,ltimes){
  
  tau <- c()
  
  for(i in 1:nrow(ltimes)){
    
    missing <- which(is.na(ltimes[i,])==TRUE)
    
    if(length(missing)!=0){
      tau[i] = sum(alpha[-missing]^2*(beta[-missing] - ltimes[i,-missing]))/sum(alpha[-missing]^2)
    } else{
      tau[i] = sum(alpha^2*(beta - ltimes[i,]))/sum(alpha^2)
    }
    
  }
  
  return(tau)
}




Lambdas <- function(ltimes, comp,alpha,beta){
  ncomp=setdiff(1: ncol(ltimes), comp)
  tcomp=PPest(alpha[comp], beta[comp], ltimes[, comp])
  tncomp=PPest(alpha[ncomp], beta[ncomp], ltimes[, ncomp])
  return((tcomp-tncomp)/sqrt(1/sum((alpha[comp])^2) + 1/sum((alpha[ncomp])^2)))
}

###############################################################################

set.seed(4102021)

TypeI <- matrix(nrow=100,ncol=4)

for(R in 1:100){

  data <- sim_dglnrt()
  
  ly        <- data.frame(log(data$rt[1:25]))

  n <- ncol(ly)

  vcov <- cov(ly,use='pairwise.complete.obs')

  model     <- paste("f1=~", 
                   paste0("a*X",1:(n-1)," + ", collapse=""), 
                   paste("a*X", n,sep=""))

  fit       <- cfa(model, 
                   sample.cov = vcov, 
                   sample.nobs = 93,
                   sample.mean = colMeans(ly,na.rm=TRUE))

  pars      <- coef(fit)


  alpha.est <- 1/sqrt(pars[1:ncol(ly)])
  beta.est  <- pars[(ncol(ly)+2):(2*ncol(ly)+1)]


  L <- Lambdas(ltimes = as.matrix(ly),
               comp   = seq(2,25,2),
               alpha  = alpha.est,
               beta   = beta.est)

  TypeI[R,] = c(length(which(L>qnorm(0.9))),
                length(which(L>qnorm(0.95))),
                length(which(L>qnorm(0.99))),
                length(which(L>qnorm(0.999))))
 
  print(R)               
}

round(colMeans(TypeI)/93,4)

round(apply(TypeI,2,min)/93,4)

round(apply(TypeI,2,max)/93,4)


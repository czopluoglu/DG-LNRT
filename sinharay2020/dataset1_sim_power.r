####################################################################
#            MODEL PARAMETER GENERATION
#
# These parameters are reported in Table 4
####################################################################

sim_dglnrt <- function() {
  
  # MODEL PARAMETERS
  
  beta  <- rnorm(25,4.19,0.38)
  alpha <- rnorm(25,1.42,0.37)
  
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
  
  C    <- c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0)
  
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

fpr <- matrix(nrow=1000,ncol=4)
tpr <- matrix(nrow=1000,ncol=4)
pre <- matrix(nrow=1000,ncol=4)
datas <- vector('list',1000)
params <- vector('list',1000)
L <- vector('list',1000)


for(R in 884:1000){
  
  data <- sim_dglnrt()
  
  datas[[R]] <- data
  
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
  
  params[[R]] <- pars
  
  alpha.est <- 1/sqrt(pars[1:ncol(ly)])
  beta.est  <- pars[(ncol(ly)+2):(2*ncol(ly)+1)]
  
  
  L[[R]] <- Lambdas(ltimes = as.matrix(ly),
               comp   = seq(2,25,2),
               alpha  = alpha.est,
               beta   = beta.est)
  
  print(R)               
}

th = .999

out <- data.frame(matrix(NA,1000,4))
colnames(out) <- c('FPR','TPR','PR_1','PR_2')

for(R in 1:1000){
  
  out[R,]$FPR  = sum(L[[R]][1:33]>qnorm(th))
  out[R,]$TPR  = sum(L[[R]][34:93]>qnorm(th))
  out[R,]$PR_1 = sum(ifelse(datas[[R]]$rt[L[[R]]>qnorm(th),]$gr==1,1,0))
  out[R,]$PR_2 = sum(ifelse(datas[[R]]$rt[L[[R]]>qnorm(th),]$gr==1,0,1))
}

# FPR across 1000 replications

sum(out$FPR)/33000  

# TPR across 1000 replications

sum(out$TPR)/60000  

# Precision across 1000 replications

sum(out$PR_2)/(sum(out$PR_1)+sum(out$PR_2))




################################################################################

# Parameter recovery

tr <- matrix(NA,1000,25)
est <- matrix(NA,1000,25)

for(R in 1:1000){
  
  tr[R,]  <- datas[[R]]$a
  est[R,] <- 1/sqrt(params[[R]][1:25])
  
}

pos1 <- seq(1,25,2)
pos2 <- seq(2,25,2)


round(mean(colMeans((est[,pos1] - tr[,pos1]))),3)
round(mean(colMeans((est[,pos2] - tr[,pos2]))),3)



tr2 <- matrix(NA,1000,25)
est2 <- matrix(NA,1000,25)

for(R in 1:1000){
  
  tr2[R,]  <- datas[[R]]$b
  est2[R,] <- params[[R]][27:51]
  
}

round(mean(colMeans((est2[,pos1] - tr2[,pos1]))),3)
round(mean(colMeans((est2[,pos2] - tr2[,pos2]))),3)




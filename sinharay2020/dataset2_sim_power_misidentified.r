require(MASS)
require(lavaan)
require(psych)

####################################################################
####################################################################
####################################################################
####################################################################
####################################################################
####################################################################

# A generic function to simulate DG-LNRT data with no item preknowledge

sim_dglnrt <- function() {
  
  # MODEL PARAMETERS
  
  beta  <- rnorm(170,3.98,0.32)
  alpha <- rnorm(170,2.06,0.28)
  
  # Tau for unflagged examinees
  
  cor0 <- matrix(c(1,.95,.95,1),2,2)
  tau0 <- mvrnorm(1590,c(0,0),cor2cov(cor0,c(0.17,0.16)))
  
  # Tau for flagged examinees 
  
  cor1 <- matrix(c(1,.92,.92,1),2,2)
  tau1 <- mvrnorm(46,c(0.19,0.28),cor2cov(cor1,c(0.35,0.30)))
  
  # A vector for item status (0: not disclosed, 1:disclosed)
  
  C    <- c(rep(0,106),rep(1,64))
  
  # Shuffle the vector, randomly selected 64 items are compromised
  
  C <- C[base::sample(1:170,170)]
  
  # RESPONSE TIME GENERATION
  
  # Note that the response time data generated is
  # already on the log scale
  
  # Unflagged Examinees
  
  rt0 <- matrix(nrow = 1590, ncol = 170)
  
  for (i in 1:1590) {
    for (j in 1:170) {
      p_t = beta[j] - tau0[i, 1]
      p_c = beta[j] - tau0[i, 2]
      p   = p_t * (1 - C[j]) + p_c * C[j]
      rt0[i, j] = rnorm(1, p, 1 / alpha[j])
    }
  }
  
  # Flagged Examinees
  
  rt1 <- matrix(nrow = 46, ncol = 170)
  
  for (i in 1:46) {
    for (j in 1:170) {
      p_t = beta[j] - tau1[i, 1]
      p_c = beta[j] - tau1[i, 2]
      p   = p_t * (1 - C[j]) + p_c * C[j]
      rt1[i, j] = rnorm(1, p, 1 / alpha[j])
    }
  }
  
  # Combine the groups
  
  rt <- rbind(cbind(data.frame(exp(rt0)), gr = 1),
              cbind(data.frame(exp(rt1)), gr = 2))
  
  
  # Random shuffle of examinees
  
  rt <- rt[sample(1:nrow(rt),nrow(rt),replace = FALSE),]
  
  return(list(
    rt = rt,
    b = beta,
    a = alpha,
    tau_t = c(tau0[, 1], tau1[, 1]),
    tau_c = c(tau0[, 2], tau1[, 2]),
    C = C
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

##############################################################################


set.seed(4102021)

datas  <- vector('list',100)
params <- vector('list',100)
L      <- vector('list',100)

for(R in 1:100){
  
  data <- sim_dglnrt()
  
  datas[[R]] <- data
  
  ly        <- data.frame(log(data$rt[,1:170]))
  
  n <- ncol(ly)
  
  
  model     <- paste("f1=~", 
                     paste0("a*X",1:(n-1)," + ", collapse=""), 
                     paste("a*X", n,sep=""))
  
  fit       <- cfa(model = model, 
                   data = ly,
                   missing = 'fiml')
  
  pars      <- coef(fit)
  
  params[[R]] <- pars
  
  alpha.est <- 1/sqrt(pars[1:ncol(ly)])
  beta.est  <- pars[(ncol(ly)+2):(2*ncol(ly)+1)]
  
  # Take a random sample of 32 compromised items out of 64
  # the remaining 32 will be treated as uncompromised
  
  truly_identified <- sample(which(data$C==1),32)
  
  # Take a random sample of 32 items out of 106 uncompromised items
  # they will be incorrectly treated as compromised
  
  misidentified <- sample(which(data$C==0),32)
  
  
  L[[R]] <- Lambdas(ltimes = as.matrix(ly),
                    comp   = c(truly_identified,misidentified),
                    alpha  = alpha.est,
                    beta   = beta.est)
  print(R)            
}

th = .999

out <- data.frame(matrix(NA,100,4))
colnames(out) <- c('FPR','TPR')

for(R in 1:100){
  
  loc1 <- which(datas[[R]]$rt$gr==1)
  loc2 <- which(datas[[R]]$rt$gr==2)
  
  out[R,]$FPR  = sum(L[[R]][loc1]>qnorm(th))
  out[R,]$TPR  = sum(L[[R]][loc2]>qnorm(th))
}

# FPR across 100 replications

sum(out$FPR)/(1590*100)  

# TPR across 1000 replications

sum(out$TPR)/(46*100)  

# Precision across 1000 replications

sum(out$TPR)/(sum(out$TPR)+sum(out$FPR))



################################################################################

# Parameter recovery

tr <- matrix(NA,100,170)
est <- matrix(NA,100,170)

for(R in 1:100){
  
  tr[R,]  <- datas[[R]]$a
  est[R,] <- 1/sqrt(params[[R]][1:170])
  
}


bias1 <- c()
bias2 <- c()

for(R in 1:100){
  
  temp <- tr[R,] - est[R,]
  
  bias1[R] = mean(temp[which(datas[[R]]$C==0)])
  
  bias2[R] = mean(temp[which(datas[[R]]$C==1)])
  
}

round(mean(bias1),3)
round(mean(bias2),3)



tr <- matrix(NA,100,170)
est <- matrix(NA,100,170)

for(R in 1:100){
  
  tr[R,]  <- datas[[R]]$b
  est[R,] <- params[[R]][255:507]
  
}


bias1 <- c()
bias2 <- c()

for(R in 1:100){
  
  temp <- tr[R,] - est[R,]
  
  bias1[R] = mean(temp[which(datas[[R]]$C==0)])
  
  bias2[R] = mean(temp[which(datas[[R]]$C==1)])
  
}

round(mean(bias1),3)
round(mean(bias2),3)
















################################################################################
################################################################################
################################################################################

# No missing data

################################################################################
################################################################################
################################################################################

# A generic function to simulate DG-LNRT data with no item preknowledge

sim_dglnrt <- function(N,pN,n,pi) {
  
  # N, total sample size
  # p, proportion of examinees with item preknowledge
  # n, number of items
  # pi,proportion of compromised items
  
  # MODEL PARAMETERS
  
  beta  <- rnorm(n,4.03,0.32)
  alpha <- rnorm(n,1.46,0.09)
  
  # Tau for unflagged examinees
  
  cor0 <- matrix(c(1,.95,.95,1),2,2)
  tau0 <- mvrnorm(N*(1-pN),c(0,0),cor2cov(cor0,c(0.17,0.16)))
  
  # Tau for flagged examinees 
  
  cor1 <- matrix(c(1,.92,.92,1),2,2)
  tau1 <- mvrnorm(N*pN,c(0.19,0.28),cor2cov(cor1,c(0.35,0.30)))
  
  # A vector for item status (0: not disclosed, 1:disclosed)
  
  C    <- rep(0,n)
  C[sample(1:n,n*pi)] = 1
  
  # RESPONSE TIME GENERATION
  
  # Note that the response time data generated is
  # already on the log scale
  
  # Unflagged Examinees
  
  rt0 <- matrix(nrow = nrow(tau0), ncol = n)
  
  for (i in 1:nrow(tau0)) {
    for (j in 1:n) {
      p_t = beta[j] - tau0[i, 1]
      p_c = beta[j] - tau0[i, 2]
      p   = p_t * (1 - C[j]) + p_c * C[j]
      rt0[i, j] = rnorm(1, p, 1 / alpha[j])
    }
  }
  
  # Flagged Examinees
  
  rt1 <- matrix(nrow = nrow(tau1), ncol = n)
  
  for (i in 1:nrow(tau1)) {
    for (j in 1:n) {
      p_t = beta[j] - tau1[i, 1]
      p_c = beta[j] - tau1[i, 2]
      p   = p_t * (1 - C[j]) + p_c * C[j]
      rt1[i, j] = rnorm(1, p, 1 / alpha[j])
    }
  }
  
  # Combine the groups
  
  rt <- rbind(cbind(data.frame(exp(rt0)), gr = 1),
              cbind(data.frame(exp(rt1)), gr = 2))
  
  
  # Random shuffle of examinees
  
  rt <- rt[sample(1:nrow(rt),nrow(rt),replace = FALSE),]
  
  return(list(
    rt = rt,
    b = beta,
    a = alpha,
    tau_t = c(tau0[, 1], tau1[, 1]),
    tau_c = c(tau0[, 2], tau1[, 2]),
    C = C
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

##############################################################################


set.seed(4102021)

fpr <- matrix(nrow=100,ncol=4)
tpr <- matrix(nrow=100,ncol=4)
pre <- matrix(nrow=100,ncol=4)
datas <- vector('list',100)
params <- vector('list',100)

for(R in 1:100){
  
  data <- sim_dglnrt(N=1650,n=170,pN=0.04,pi=0.3)
  
  datas[[R]] <- data
  
  ly        <- data.frame(log(data$rt[,1:(ncol(data$rt)-1)]))
  
  n <- ncol(ly)
  
  
  model     <- paste("f1=~", 
                     paste0("a*X",1:(n-1)," + ", collapse=""), 
                     paste("a*X", n,sep=""))
  
  fit       <- cfa(model = model, 
                   data = ly,
                   meanstructure=TRUE,
                   auto.var=TRUE)
  
  pars      <- coef(fit)
  params[[R]] <- pars
  
  alpha.est <- 1/sqrt(pars[1:ncol(ly)])
  beta.est  <- pars[(ncol(ly)+2):(2*ncol(ly)+1)]
  
  
  L <- Lambdas(ltimes = as.matrix(ly),
               comp   = which(data$C==1),
               alpha  = alpha.est,
               beta   = beta.est)
  
  fpr[R,] = c(sum(L[which(data$rt$gr==1)]>qnorm(0.9))/length(which(data$rt$gr==1)),
              sum(L[which(data$rt$gr==1)]>qnorm(0.95))/length(which(data$rt$gr==1)),
              sum(L[which(data$rt$gr==1)]>qnorm(0.99))/length(which(data$rt$gr==1)),
              sum(L[which(data$rt$gr==1)]>qnorm(0.999))/length(which(data$rt$gr==1)))
  
  tpr[R,] = c(sum(L[which(data$rt$gr==2)]>qnorm(0.9))/length(which(data$rt$gr==2)),
              sum(L[which(data$rt$gr==2)]>qnorm(0.95))/length(which(data$rt$gr==2)),
              sum(L[which(data$rt$gr==2)]>qnorm(0.99))/length(which(data$rt$gr==2)),
              sum(L[which(data$rt$gr==2)]>qnorm(0.999))/length(which(data$rt$gr==2)))
  
  pre[R,] = c(mean(ifelse(data$rt[which(L>qnorm(0.9)),]$gr==1,0,1)),
              mean(ifelse(data$rt[which(L>qnorm(0.95)),]$gr==1,0,1)),
              mean(ifelse(data$rt[which(L>qnorm(0.99)),]$gr==1,0,1)),
              mean(ifelse(data$rt[which(L>qnorm(0.999)),]$gr==1,0,1)))
  
  
  print(R)            
}


round(colMeans(fpr),3)
round(apply(fpr,2,min),3)
round(apply(fpr,2,max),3)

round(colMeans(tpr),3)
round(apply(tpr,2,min),3)
round(apply(tpr,2,max),3)

round(colMeans(pre),3)
round(apply(pre,2,min),3)
round(apply(pre,2,max),3)


################################################################################

# Parameter recovery

tr <- matrix(NA,100,170)
est <- matrix(NA,100,170)

for(R in 1:100){
  
  tr[R,]  <- datas[[R]]$a
  est[R,] <- 1/sqrt(params[[R]][1:170])
  
}


bias1 <- c()
bias2 <- c()

for(R in 1:100){
  
  temp <- tr[R,] - est[R,]
  
  bias1[R] = mean(temp[which(datas[[R]]$C==0)])
  
  bias2[R] = mean(temp[which(datas[[R]]$C==1)])
  
}

round(mean(bias1),3)
round(mean(bias2),3)



tr <- matrix(NA,100,170)
est <- matrix(NA,100,170)

for(R in 1:100){
  
  tr[R,]  <- datas[[R]]$b
  est[R,] <- params[[R]][172:341]
  
}


bias1 <- c()
bias2 <- c()

for(R in 1:100){
  
  temp <- tr[R,] - est[R,]
  
  bias1[R] = mean(temp[which(datas[[R]]$C==0)])
  
  bias2[R] = mean(temp[which(datas[[R]]$C==1)])
  
}

round(mean(bias1),3)
round(mean(bias2),3)



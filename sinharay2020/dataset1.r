
load(here("data/dglnrt_v1/results with original priors.RData"))

head(d.long)
head(d.sub)

rt <- log(d.sub[,3:27])
    
ly        <- data.frame(rt)

n <- ncol(ly)

vcov <- cov(ly,use='pairwise.complete.obs')
  
  
  
model     <- paste("f1=~", 
                   paste0("a*Q",1:(n-1),"RT + ", collapse=""), 
                   paste("a*Q", n,"RT",sep=""))

fit       <- cfa(model = model, 
                 data = ly,
                 missing = 'fiml')

pars      <- coef(fit)


alpha.est <- 1/sqrt(pars[1:ncol(ly)])
beta.est  <- pars[(ncol(ly)+2):(2*ncol(ly)+1)]



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
  tall=PPest(alpha, beta,ltimes)
  return((tcomp-tncomp)/sqrt(1/sum((alpha[comp])^2) + 1/sum((alpha[ncomp])^2)))
}


L <- Lambdas(ltimes = as.matrix(rt),
             comp   = seq(2,25,2),
             alpha  = alpha.est,
             beta   = beta.est)

out <- cbind(L,ifelse(d.sub$COND==1,0,1))

out  <- na.omit(out)

length(which(out[which(out[,2]==0),1]>qnorm(0.9)))/30
length(which(out[which(out[,2]==0),1]>qnorm(0.95)))/30
length(which(out[which(out[,2]==0),1]>qnorm(0.99)))/30
length(which(out[which(out[,2]==0),1]>qnorm(0.999)))/30


length(which(out[which(out[,2]==1),1]>qnorm(0.9)))/60
length(which(out[which(out[,2]==1),1]>qnorm(0.95)))/60
length(which(out[which(out[,2]==1),1]>qnorm(0.99)))/60
length(which(out[which(out[,2]==1),1]>qnorm(0.999)))/60


table(out[which(out[,1]>qnorm(0.9)),2])
table(out[which(out[,1]>qnorm(0.95)),2])
table(out[which(out[,1]>qnorm(0.99)),2])
table(out[which(out[,1]>qnorm(0.999)),2])


load(here("data/dglnrt_v2/results.RData"))

d.wide <- reshape(d.long,
                  idvar = 'id',
                  v.names = 'RT',
                  timevar = 'item',
                  direction = 'wide')

table(d.wide$p_flag)


i.st <- c()

for(i in 1:253){
  
  i.st[i] = unique(d.long[which(d.long$item==i),]$i.status)
  
}

comp <- which(i.st==1)

################################################################################

rt <- d.wide[,8:260]

ilabels <- as.numeric(substring(colnames(rt),4))

rt <- rt[,order(ilabels)]

ly        <- data.frame(rt)

n <- ncol(ly)

model     <- paste("f1=~", 
                   paste0("a*RT.",1:(n-1)," + ", collapse=""), 
                   paste("a*RT.", n,sep=""))

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
  return((tcomp-tncomp)/sqrt(1/sum((alpha[comp])^2) + 1/sum((alpha[ncomp])^2)))
}


L <- Lambdas(ltimes = as.matrix(rt),
             comp   = comp,
             alpha  = alpha.est,
             beta   = beta.est)

out <- cbind(L,d.wide$p_flag)

out  <- na.omit(out)

length(which(out[which(out[,2]==0),1]>qnorm(0.9)))/3186
length(which(out[which(out[,2]==0),1]>qnorm(0.95)))/3186
length(which(out[which(out[,2]==0),1]>qnorm(0.99)))/3186
length(which(out[which(out[,2]==0),1]>qnorm(0.999)))/3186


length(which(out[which(out[,2]==1),1]>qnorm(0.9)))/94
length(which(out[which(out[,2]==1),1]>qnorm(0.95)))/94
length(which(out[which(out[,2]==1),1]>qnorm(0.99)))/94
length(which(out[which(out[,2]==1),1]>qnorm(0.999)))/94


table(out[which(out[,1]>qnorm(0.9)),2])
table(out[which(out[,1]>qnorm(0.95)),2])
table(out[which(out[,1]>qnorm(0.99)),2])
table(out[which(out[,1]>qnorm(0.999)),2])























require(here)
require(rstan)
require(psych)
require(ggplot2)
require(gridExtra)
require(qqplotr)
require(reshape2)

################################################################################
# Load the pre-saved model output

load(here("data/dglnrt_v2/results.RData"))

################################################################################
d.wide <- reshape(d.long[,c('id','RT','item','p_flag')],
                  idvar = 'id',
                  v.names = 'RT',
                  timevar = 'item',
                  direction = 'wide')

table(d.wide$p_flag)


i.st <- c()

for(i in 1:253){
  
  i.st[i] = unique(d.long[which(d.long$item==i),]$i.status)
  
}

comp  <- which(i.st==1)
comp2 <- which(i.st==0)

p.st <- c()

for(i in 1:3280){
  
  p.st[i] = unique(d.long[which(d.long$id==i),]$p_flag)
  
}

cheaters  <- which(p.st==1)
cheaters2 <- which(p.st==0)



rt <- d.wide[,3:255]

ilabels <- as.numeric(substring(colnames(rt),4))

rt <- rt[,order(ilabels)]


################################################################################
# Extract the betas

betas <- summary(stanfit, pars = c("beta"), probs = c(0.025, 0.975))$summary
betas 
describe(betas[,1])

# Extract the alphas

alphas <- summary(stanfit, pars = c("alpha"), probs = c(0.025, 0.975))$summary
alphas
describe(alphas[,1]^2)


# Extract the tau_t

tau_t <- summary(stanfit, pars = c("tau_t"), probs = c(0.025, 0.975))$summary
tau_t
describe(tau_t[,1])

# Extract the tau_c

tau_c <- summary(stanfit, pars = c("tau_c"), probs = c(0.025, 0.975))$summary
tau_c
describe(tau_c[,1])

# Extract Ts

Ts <- summary(stanfit, pars = c("T"), probs = c(0.025, 0.975))$summary
Ts
################################################################################

# Correlation between Observed vs Predicted

betas <- betas[,1]
alphas <- alphas[,1]^2
tau_t <- tau_t[,1]
tau_c <- tau_c[,1]
Ts <- Ts[,1]


rt_pred <- rt


rt_pred[cheaters,comp] = 
  matrix(betas[comp],nrow = length(cheaters),ncol=length(comp),byrow=TRUE) -
  matrix(tau_c[cheaters],nrow=length(cheaters),ncol=length(comp),byrow=FALSE) 

rt_pred[cheaters,comp2] = 
  matrix(betas[comp2],nrow = length(cheaters),ncol=length(comp2),byrow=TRUE) -
  matrix(tau_t[cheaters],nrow=length(cheaters),ncol=length(comp2),byrow=FALSE)

rt_pred[cheaters2,comp] = 
  matrix(betas[comp],nrow = length(cheaters2),ncol=length(comp),byrow=TRUE) -
  matrix(tau_t[cheaters2],nrow=length(cheaters2),ncol=length(comp),byrow=FALSE)

rt_pred[cheaters2,comp2] = 
  matrix(betas[comp2],nrow = length(cheaters2),ncol=length(comp2),byrow=TRUE) -
  matrix(tau_t[cheaters2],nrow=length(cheaters2),ncol=length(comp2),byrow=FALSE)


rt.l      <- stack(rt)
rt_pred.l <- stack(rt_pred)


cor(rt.l$values, rt_pred.l$values,use='pairwise')

loc <- sample(1:nrow(rt.l),10000)

p <- ggplot()+
  geom_point(aes(x=rt.l[loc,]$values, y = rt_pred.l[loc,]$values))+
  theme_bw()+
  ggtitle('Correlation between Observed Log Response Time and Predicted Log Response Time')+
  xlab('Observed Log Response Time')+
  ylab('Predicted Log Response Time')
  
  
  

ggsave(filename = 'Correlation.jpeg',
       plot     = p,
       device   = 'jpeg',
       path     = here('modelfit/dataset2'),
       width    = 13,
       height   = 8)


################################################################################

# Scatteer Plot, Latent Speed vs Standardized Residual
# Compromised ITems

pl <- lapply (comp,
              function(item){
                
                temp     <- as.data.frame(rt[,item] - rt_pred[,item])
                temp$tau <- tau_t
                temp$res <- temp[,1]/(1/alphas[item])
                
                ggplot(data=temp,aes(x=tau,y=res))+
                  geom_point()+
                  theme_bw() + 
                  xlab('Latent Speed Parameter')+
                  ylab('Standardized Residual')+
                  ggtitle(paste0("Item ",item))+
                  ylim(c(-5,5))+
                  geom_hline(yintercept=-2,lty=2)+
                  geom_hline(yintercept=2,lty=2)
              }
)


ind <- seq(1,length(pl)+1,6)

for(i in 1:(length(ind)-1)){
  
  ggsave(filename = paste0('fit',i,'.jpeg'),
       plot     = grid.arrange(grobs=pl[ind[i]:(ind[i+1]-1)], nrow=2, ncol=3),
       device   = 'jpeg',
       path     = here('modelfit/dataset2'),
       width    = 13,
       height   = 8)
  
  print(i)
}



# Uncompromised ITems

pl <- lapply (comp2,
              function(item){
                
                temp     <- as.data.frame(rt[,item] - rt_pred[,item])
                temp$tau <- tau_t
                temp$res <- temp[,1]/(1/alphas[item])
                
                ggplot(data=temp,aes(x=tau,y=res))+
                  geom_point()+
                  theme_bw() + 
                  xlab('Latent Speed Parameter')+
                  ylab('Standardized Residual')+
                  ggtitle(paste0("Item ",item))+
                  ylim(c(-5,5))+
                  geom_hline(yintercept=-2,lty=2)+
                  geom_hline(yintercept=2,lty=2)
              }
)


ind <- seq(1,length(pl)+1,6)

for(i in 1:(length(ind)-1)){
  
  ggsave(filename = paste0('fit',i+15,'.jpeg'),
         plot     = grid.arrange(grobs=pl[ind[i]:(ind[i+1]-1)], nrow=2, ncol=3),
         device   = 'jpeg',
         path     = here('modelfit/dataset2'),
         width    = 13,
         height   = 8)
  
  print(i)
}

################################################################################

res.wide <- rt - rt_pred

for(i in 1:253){
  
  res.wide[,i]/(1/alphas[i])
  
}


boxplot(res.wide,col="gray",xaxt="n",ylab="Standardized Residual",xlab="Item Location")
axis(side=1,at=1:253)

describe(colMeans(res.wide[,comp],na.rm=TRUE))
hist(colMeans(res.wide[,comp],na.rm=TRUE))


describe(colMeans(res.wide[,comp2],na.rm=TRUE))
hist(colMeans(res.wide[,comp2],na.rm=TRUE))


plot(1:253,colMeans(res.wide,na.rm=TRUE),ylim=c(-.2,.2),type="l",
     ylab="Mean Standardized Residual",xlab="Item Location",xaxt="n")


abline(h=0,lty=2,col="gray")

res.cor <- cor(res.wide,use="pairwise.complete.obs")


diag(res.cor) <- NA
describe(res.cor)




melted_cormat <- melt(res.cor)

ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()



# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(res.cor)
upper_tri
melted_cormat <- melt(upper_tri, na.rm=TRUE)

dim(melted_cormat)

melted_cormat[,1] = as.numeric(substring(melted_cormat[,1],first=4))
melted_cormat[,2] = as.numeric(substring(melted_cormat[,2],first=4))

melted_cormat = melted_cormat[which(abs(melted_cormat[,3])>.12),]


dim(melted_cormat)

p <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2)) + 
  geom_tile(color="black")+theme_bw()+xlab("Item Location")+ylab("Item Location")+
  scale_fill_continuous(name = "")+
  ggtitle('Residual Correlations (>0.12)')

ggsave(filename = 'Residual Correlations.jpeg',
       plot     = p,
       device   = 'jpeg',
       path     = here('modelfit/dataset2'),
       width    = 13,
       height   = 8)



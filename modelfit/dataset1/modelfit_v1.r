
require(here)
require(rstan)
require(psych)
require(ggplot2)
require(gridExtra)
require(qqplotr)
require(reshape2)

################################################################################

# Load the pre-saved model output

load(here('data/dglnrt_v1/results with original priors.RData'))

d2 <- read.csv(here('data/dglnrt_v1/UVA Experiment Data.csv'),na.strings = '#N/A')

################################################################################

pl <- lapply (5:29,
              function(item){

                
                temp <- d[,c(1:4,item)]          
                colnames(temp)[5] <- 'y'
                temp$y <- log(temp$y)
                
                ggplot(data=temp,aes(sample=y))+
                  stat_qq_band(alpha=.2)+
                  stat_qq_line()+
                  stat_qq_point()+
                  theme_bw()+
                  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
                  ggtitle(paste0("Item ",item-4))
                  
              }
)

grid.arrange(grobs=pl[1:6], nrow=2, ncol=3)
grid.arrange(grobs=pl[7:12], nrow=2, ncol=3)
grid.arrange(grobs=pl[13:18], nrow=2, ncol=3)
grid.arrange(grobs=pl[19:24], nrow=2, ncol=3)
grid.arrange(grobs=pl[25], nrow=2, ncol=3)
################################################################################

# Extract the betas

betas <- summary(stanfit, pars = c("beta"), probs = c(0.025, 0.975))$summary
betas 
describe(betas[,1])

# Extract the alphas

alphas <- summary(stanfit, pars = c("alpha"), probs = c(0.025, 0.975))$summary
alphas
describe(alphas[,1])


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
alphas <- alphas[,1]
tau_t <- tau_t[,1]
tau_c <- tau_c[,1]
Ts <- Ts[,1]

d.long$pred <- NA

for (i in 1:nrow(d.long)) {
  
  pt <- betas[d.long[i,]$Item] - tau_t[d.long[i,]$ID]
  pc <- betas[d.long[i,]$Item] - tau_c[d.long[i,]$ID]
  
  p = (pt^(1-Ts[d.long[i,]$ID]))*
      (((1-d.long[i,]$i.status)*pt + 
          (d.long[i,]$i.status)*pc)^Ts[d.long[i,]$ID])
  
  d.long[i,]$pred <- p
}

cor(d.long$logRT,d.long$pred)

plot(d.long$logRT,d.long$pred,
     xlab='Observed Log Response Time',
     ylab='Predicted Log Response Time',)


################################################################################

# Scatteer Plot, Latent Speed vs Standardized Residual

pl <- lapply (seq(1,25,2),
              function(item){
                
                temp <- d.long[which(d.long$Item==item),]
                temp$tau <- tau_t[temp$ID]
                temp$res <- (temp$logRT - temp$pred)/(1/alphas[item])
                                         
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

grid.arrange(grobs=pl[1:6], nrow=2, ncol=3)
grid.arrange(grobs=pl[7:13], nrow=2, ncol=4)




pl <- lapply (seq(2,25,2),
              function(item){
                
                temp <- d.long[which(d.long$Item==item),]
                tau  <- ifelse(temp$COND==1,tau_t,tau_c)
                temp$tau <- tau[temp$ID]
                temp$res <- (temp$logRT - temp$pred)/(1/alphas[item])
                
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

grid.arrange(grobs=pl[1:6], nrow=2, ncol=3)
grid.arrange(grobs=pl[7:12], nrow=2, ncol=3)

################################################################################


d.long$sres <- NA

for(i in 1:25){
  
  loc <- which(d.long$Item==i)
  
  d.long[loc,]$sres = (d.long[loc,]$logRT - d.long[loc,]$pred)/(1/alphas[i])
  
}


res.wide <- reshape(data=d.long[,c('ID','Item','sres')],
                    idvar="ID",
                    v.names="sres",
                    timevar="Item",
                    direction="wide")

res.wide <- res.wide[,-1]

boxplot(res.wide,col="gray",xaxt="n",ylab="Standardized Residual",xlab="Item Location")
axis(side=1,at=1:25)

plot(1:25,colMeans(res.wide,na.rm=TRUE),ylim=c(-.2,.2),type="l",
     ylab="Mean Standardized Residual",xlab="Item Location",xaxt="n")

axis(side=1,at=1:25)

abline(h=0,lty=2,col="gray")

res.cor <- cor(res.wide,use="pairwise.complete.obs")

res.cor[1,2]
res.cor[2,3]
res.cor[3,4]
res.cor[4,5]
res.cor[5,6]
res.cor[6,7]
res.cor[7,8]
res.cor[8,9]
res.cor[9,10]
res.cor[10,11]
res.cor[11,12]
res.cor[12,13]
res.cor[13,14]
res.cor[14,15]
res.cor[15,16]
res.cor[16,17]
res.cor[17,18]
res.cor[18,19]
res.cor[19,20]
res.cor[20,21]
res.cor[21,22]
res.cor[22,23]
res.cor[23,24]
res.cor[24,25]

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

melted_cormat[,1] = as.numeric(substring(melted_cormat[,1],first=6))
melted_cormat[,2] = as.numeric(substring(melted_cormat[,2],first=6))

melted_cormat = melted_cormat[which(abs(melted_cormat[,3])>.38),]

ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(color="black")+theme_bw()+xlab("Item Location")+ylab("Item Location")+
  scale_x_continuous(breaks=1:25,limits=c(.5,25.5))+
  scale_y_continuous(breaks=1:25,limits=c(.5,25.5))+
  scale_fill_continuous(name = "")+
  scale_fill_gradientn(colours = c("white", "black"), values = c(0,1))










































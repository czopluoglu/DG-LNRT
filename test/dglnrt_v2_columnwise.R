require(cmdstanr)
require(here)
require(psych)
########################

# This fits the model to the real dataset by Cizek and Wollack

#########################


d.long <- read.csv(here('data/rt_long.csv'))

d.long <- na.omit(d.long)

d.long$item <- as.numeric(substring(d.long$Item,6))

icode = c(201,204,205,206,208,210,211,214,216,218,220,222,223,224,226,227,228,
          229,230,231,232,234,236,237,239,241,242,243,244,247,248,249,253,255,
          257,258,261,263,264,265,267,268,269,272,279,280,281,282,283,286,287,
          288,292,296,297,298,302,303,305,307,308,314,316,318,322,323,327,329,
          331,338,340,341,343,345,346,347,348,350,353,359,366,367,368)

icode2 = 171:253


for(i in 1:83){
  d.long[which(d.long$item == icode[i]),]$item = icode2[i]
}

unique(d.long$item)

tab <- as.data.frame(table(d.long$Item,d.long$item))
tab = tab[which(tab[,3]!=0),]



d.long$id <- NA

ind = which(substring(d.long$EID,1,3)=='e10')
d.long[ind,]$id = as.numeric(substring(d.long[ind,]$EID,4,7))

ind = which(substring(d.long$EID,1,3)=='e20')
d.long[ind,]$id = as.numeric(substring(d.long[ind,]$EID,4,7))+1636

tab <- as.data.frame(table(d.long$EID,d.long$id))
tab = tab[which(tab[,3]!=0),]


unique(d.long$id)

flag.key <- c()
for(i in 1:253){
  flag.key[i] = unique(d.long[which(d.long$item==i),]$i.status)
}

n.key <- c()
for(i in 1:253){
  n.key[i] = length(which(d.long$item==i))
}
  
########################################################################

d.long$i.status <- ifelse(d.long$i_flag=='Flagged',1,0)


d.long <- d.long[order(d.long$item),]
#########################################################################


p1 <- d.long[which(d.long$item==1),]$id
y1 <- d.long[which(d.long$item==1),]$RT
i1 <- d.long[which(d.long$item==1),]$i.status
n1 <- length(y1) 
  
p2 <- d.long[which(d.long$item==2),]$id
y2 <- d.long[which(d.long$item==2),]$RT
i2 <- d.long[which(d.long$item==2),]$i.status
n2 <- length(y1)

p3 <- d.long[which(d.long$item==3),]$id
y3 <- d.long[which(d.long$item==3),]$RT
i3 <- d.long[which(d.long$item==3),]$i.status
n3 <- length(y3)

p4 <- d.long[which(d.long$item==4),]$id
y4 <- d.long[which(d.long$item==4),]$RT
i4 <- d.long[which(d.long$item==4),]$i.status
n4 <- length(y4)

p5 <- d.long[which(d.long$item==5),]$id
y5 <- d.long[which(d.long$item==5),]$RT
i5 <- d.long[which(d.long$item==5),]$i.status
n5 <- length(y5)

p6 <- d.long[which(d.long$item==6),]$id
y6 <- d.long[which(d.long$item==6),]$RT
i6 <- d.long[which(d.long$item==6),]$i.status
n6 <- length(y6)

p7 <- d.long[which(d.long$item==7),]$id
y7 <- d.long[which(d.long$item==7),]$RT
i7 <- d.long[which(d.long$item==7),]$i.status
n7 <- length(y7)

p8 <- d.long[which(d.long$item==8),]$id
y8 <- d.long[which(d.long$item==8),]$RT
i8 <- d.long[which(d.long$item==8),]$i.status
n8 <- length(y8)

p9 <- d.long[which(d.long$item==9),]$id
y9 <- d.long[which(d.long$item==9),]$RT
i9 <- d.long[which(d.long$item==9),]$i.status
n9 <- length(y9)

p10 <- d.long[which(d.long$item==10),]$id
y10 <- d.long[which(d.long$item==10),]$RT
i10 <- d.long[which(d.long$item==10),]$i.status
n10 <- length(y10)

p11 <- d.long[which(d.long$item==21),]$id
y11 <- d.long[which(d.long$item==21),]$RT
i11 <- d.long[which(d.long$item==21),]$i.status
n11 <- length(y11)

p12 <- d.long[which(d.long$item==22),]$id
y12 <- d.long[which(d.long$item==22),]$RT
i12 <- d.long[which(d.long$item==22),]$i.status
n12 <- length(y12)

p13 <- d.long[which(d.long$item==30),]$id
y13 <- d.long[which(d.long$item==30),]$RT
i13 <- d.long[which(d.long$item==30),]$i.status
n13 <- length(y13)

p14 <- d.long[which(d.long$item==32),]$id
y14 <- d.long[which(d.long$item==32),]$RT
i14 <- d.long[which(d.long$item==32),]$i.status
n14 <- length(y14)

p15 <- d.long[which(d.long$item==34),]$id
y15 <- d.long[which(d.long$item==34),]$RT
i15 <- d.long[which(d.long$item==34),]$i.status
n15 <- length(y15)



#########################################################################

data_rt <- list(I              = 20,
                J              = 3280,
                n1 = n1,
                n2 = n2,
                n3 = n3,
                n4 = n4,
                n5 = n5,
                n6 = n6,
                n7 = n7,
                n8 = n8,
                n9 = n9,
                n10 = n10,
                n11 = n11,
                n12 = n12,
                n13 = n13,
                n14 = n14,
                n15 = n15,
                y1 = y1,
                y2 = y2,
                y3 = y3,
                y4 = y4,
                y5 = y5,
                y6 = y6,
                y7 = y7,
                y8 = y8,
                y9 = y9,
                y10 = y10,
                y11 = y11,
                y12 = y12,
                y13 = y13,
                y14 = y14,
                y15 = y15,
                p1 = p1,
                p2 = p2,
                p3 = p3,
                p4 = p4,
                p5 = p5,
                p6 = p6,
                p7 = p7,
                p8 = p8,
                p9 = p9,
                p10 = p10,
                p11 = p11,
                p12 = p12,
                p13 = p13,
                p14 = p14,
                p15 = p15)

########################################################################

mod <- cmdstan_model(here('dglnrt_v2_columnwise.stan'))

fit <- mod$sample(
  data = data_rt,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_warmup   = 100,
  iter_sampling = 100,
  refresh = 10,
  adapt_delta = 0.99,
  max_treedepth = 15
)

fit$cmdstan_summary()

stanfit <- rstan::read_stan_csv(fit$output_files())


summary(stanfit, pars = c("mu1","sigma1","sigma_t","sigma_c"), 
        probs = c(0.025, 0.975))$summary

betas <- summary(stanfit, pars = c("beta"), probs = c(0.025, 0.975))$summary
betas 
describe(betas[,1])

alphas <- summary(stanfit, pars = c("alpha"), probs = c(0.025, 0.975))$summary
alphas
describe(alphas[,1])

tau_t <- summary(stanfit, pars = c("tau_t"), probs = c(0.025, 0.975))$summary
tau_t
describe(tau_t[,1])

tau_c <- summary(stanfit, pars = c("tau_c"), probs = c(0.025, 0.975))$summary
tau_c
describe(tau_c[,1])

T <- summary(stanfit, pars = c("T"), probs = c(0.025, 0.975))$summary
tau_c
describe(tau_c[,1])


as.numeric(which(T[,1]>.95))
unique(d.long[which(d.long$p_flag==1),]$id)

table(unique(d.long[which(d.long$p_flag==1),]$id) %in% as.numeric(which(T[,1]>.95)))



length(unique(d.long$id))
length(unique(d.long[which(d.long$p_flag==1),]$id))












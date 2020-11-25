require(MASS)
require(MBESS)
require(matrixStats)
require(cmdstanr)
require(here)
require(rstan)
require(psych)
###############################################################################

# This is a code being used for simulating a null data (where there is no item 
# preknowledge). There are 3280 examinees and 253 items as in Real Dataset 2.
# DG-LNRT model is fitted as usual. We expect DG-LNRT to detect nobody in this 
# dataset.

# Parameters come from EMIP paper where a multigroup lognormal response time
# model was fitted with gated mechanism using Real Dataset 2.

# For every replication, we generate a new set of model parameters from the 
# the defined distributions.

# Prepare an Appendix describing the model, model estimation, fitting, etc.
# Or, put the whole paper as an appendix. 

####################################################################
####################################################################
####################################################################

# A generic function to simulate DG-LNRT data with no item preknowledge

  sim_dglnrt <- function() {
    
    # MODEL PARAMETERS
    
    beta  <- rnorm(253,3.98,0.32)
    alpha <- rnorm(253,2.06,0.28)
    
    # Tau for unflagged examinees
    
    cor0 <- matrix(c(1,.95,.95,1),2,2)
    tau0 <- mvrnorm(3186,c(0,0),cor2cov(cor0,c(0.17,0.16)))
    
    # Tau for flagged examinees 
    
    cor1 <- matrix(c(1,.92,.92,1),2,2)
    tau1 <- mvrnorm(94,c(0.19,0.28),cor2cov(cor1,c(0.35,0.30)))
       
    # A vector for item status (0: not disclosed, 1:disclosed)
    
    C    <- rep(0,253)
	
	# null condition, all zeros
    
    # RESPONSE TIME GENERATION
    
    # Note that the response time data generated is
    # already on the log scale
    
    # Unflagged Examinees
    
    rt0 <- matrix(nrow = 3186, ncol = 253)
    
    for (i in 1:3186) {
      for (j in 1:253) {
        p_t = beta[j] - tau0[i, 1]
        p_c = beta[j] - tau0[i, 2]
        p   = p_t * (1 - C[j]) + p_c * C[j]
        rt0[i, j] = rnorm(1, p, 1 / alpha[j])
      }
    }
    
    # Flagged Examinees
    
    rt1 <- matrix(nrow = 94, ncol = 253)
    
    for (i in 1:94) {
      for (j in 1:253) {
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
    
    # Introduce missingness as in the dataset 
    # Each half of the examinees respond only to 170 items
    # 87 items are in common.
    
      rt[1:1640,171:253] = NA
      rt[1641:3280,88:170] = NA
    
    
    return(list(
      rt = rt,
      b = beta,
      a = alpha,
      tau_t = c(tau0[, 1], tau1[, 1]),
      tau_c = c(tau0[, 2], tau1[, 2])
    ))
  }

##############################################################################

mod <- cmdstan_model('/gpfs/projects/edquant/cengiz/dglnrt_null/dglnrt2.stan')

    data <- sim_dglnrt()
    
    d.sub <- data$rt[,1:253]
    
    d.sub$ID <- 1:3280
    
    d.long <- reshape(
      data        = d.sub,
      idvar       = "ID",
      varying     = list(colnames(d.sub)[1:253]),
      timevar     = "Item",
      times       = 1:253,
      v.names      = "RT",
      direction   = "long"
    )
    
    d.long <- na.omit(d.long)
    d.long$logRT <- log(d.long$RT)
    
    d.long$i.status <- NA
    
    flagged.item <- sample(1:253,90)
    
    d.long[which(d.long$Item%in%flagged.item==TRUE),]$i.status = 1
    d.long[which(d.long$Item%in%flagged.item==FALSE),]$i.status = 0
    
    data_rt <- list(
      J              = 253,
      I              = 3280,
      n_obs          = length(d.long$RT),
      ind_person_obs = d.long$ID,
      ind_item_obs   = d.long$Item,
      i_status_obs   = d.long$i.status,
      Y              = d.long$logRT
    )
    
    
    fit <- mod$sample(
      data = data_rt,
      seed = 1234,
      chains = 4,
      parallel_chains = 4,
      iter_warmup   = 1000,
      iter_sampling = 2000,
      refresh = 1000,
      adapt_delta = 0.99)

    
    fit$cmdstan_summary()
    
    stanfit[[i]] <- rstan::read_stan_csv(fit$output_files())
    
    
save.image('/gpfs/projects/edquant/cengiz/dglnrt_null/rep54.RData')

# SPECIFY THE JAGS MODEL TO FIT TO YUKON CHINOOK DATA

# file to store the JAGS model code
jags_model_file = "ssm-model.txt"

# JAGS model code
jags_model = function() {
  
  # priors: population-specific parameters
  # and quantities derived from them
  for (s in 1:ns) {
    phi[s] ~ dunif(-0.99,0.99)
    U_msy[s] ~ dunif(0.01,0.99)
    log_S_msy[s] ~ dnorm(0,0.001) %_% I(1,11.5)
    S_msy[s] <- exp(log_S_msy[s])
    alpha[s] <- exp(U_msy[s])/(1 - U_msy[s])
    log_alpha[s] <- log(alpha[s])
    beta[s] <- U_msy[s]/S_msy[s]
    log_resid_0[s] <- 0
  }
  
  # priors: variance-covariance matrix on recruitment white noise
  Tau_R[1:ns,1:ns] ~ dwish(R_wish[1:ns,1:ns],df_wish)
  Sigma_R[1:ns,1:ns] <- inverse(Tau_R)
  
  # each population's recruitment SD
  for (s in 1:ns) {
    sigma_R[s] <- sqrt(Sigma_R[s,s])
  }
  
  # correlation matrix
  for (i in 1:ns) {
    for (j in 1:ns) {
      rho_mat[i,j] <- Sigma_R[i,j]/(sigma_R[i] * sigma_R[j])
    }
  }
  
  # generate expected values for recruitment states
  for (s in 1:ns) {
    # initialize population with expected recruitment at the unfished level
    R_eq[s] <- log_alpha[s]/beta[s]
    R0[s] <- R_eq[s]
    log_R0[s] <- log(R0[s])
    log_R_mean1[1,s] <- log_R0[s]
    R_mean1[1,s] <- R0[s]
    log_R_mean2[1,s] <- log_R_mean1[1,s] + phi[s] * log_resid_0[s]
    
    # years without spawner linkage: expected value is unfished recruitment
    # recruits generated by these early years become the spawners in the early years of the next section
    for (y in 2:a_max) {
      R_mean1[y,s] <- R0[s]
      log_R_mean1[y,s] <- log_R0[s]
      log_R_mean2[y,s] <- log_R_mean1[y,s] + phi[s] * log_resid[y - 1,s]
    }
    
    # years with spawner linkage: expected value is from Ricker dynamics
    for (y in (a_max + 1):ny) {
      R_mean1[y,s] <- S[y - a_max,s] * exp(log_alpha[s] - beta[s] * S[y - a_max,s])
      log_R_mean1[y,s] <- log(R_mean1[y,s])
      log_R_mean2[y,s] <- log_R_mean1[y,s] + phi[s] * log_resid[y - 1,s]
    }
  }
  
  # draw stochastic recruitment states
  # one year across populations is a multivariate normal random vector
  for (y in 1:ny) {
    log_R[y,1:ns] ~ dmnorm(log_R_mean2[y,1:ns],Tau_R[1:ns,1:ns])
  }
  
  # calculate recruitment residuals by population and year
  for (y in 1:ny) {
    for (s in 1:ns) {
      R[y,s] <- exp(log_R[y,s])
      log_resid[y,s] <- log_R[y,s] - log_R_mean1[y,s]
    }
  }
  
  # maturity schedule: assumed to be shared exactly among all populations
  # but have Dirichlet-distributed maturity vectors for each brood year
  
  # expected probabilities: mean of of hyperdistribution
  prob[1] ~ dbeta(1,1)
  prob[2] ~ dbeta(1,1)
  prob[3] ~ dbeta(1,1)
  pi[1] <- prob[1]
  pi[2] <- prob[2] * (1 - pi[1])
  pi[3] <- prob[3] * (1 - pi[1] - pi[2])
  pi[4] <- 1 - pi[1] - pi[2] - pi[3]
  
  # (inverse) dispersion parameter
  D_scale ~ dunif(0.03,1)
  D_sum <- 1/D_scale^2
  
  # draw dirichlet random vectors as a mixture of gammas
  # see Fleischman et al. (2013; 10.1139/cjfas-2012-0112)
  for (a in 1:na) {
    dir_alpha[a] <- D_sum * pi[a]
    for (y in 1:ny) {
      g[y,a] ~ dgamma(dir_alpha[a],1)
      p[y,a] <- g[y,a]/sum(g[y,1:na])
    }
  }
  
  # apportion brood year recruitment to calendar year return by age
  # for each population
  for (s in 1:ns) {
    for (t in 1:nt) {
      for (a in 1:na) {
        # abundance returning in calendar (observation) year t at the a-th possible age
        N_tas[t,a,s] <- R[t+na-a,s] * p[t+na-a,a]
        
        # proportion of the run returning to this population that is of each age
        q_stock[t,a,s] <- N_tas[t,a,s]/N[t,s]
      }
      
      # calculate totals across age for each population
      N[t,s] <- sum(N_tas[t,1:na,s])
      S[t,s] <- N[t,s] * (1 - U[t] * v[s])
      C[t,s] <- N[t,s] * (U[t] * v[s])
    }
  }
  
  # calculate totals across populations including age composition
  for (t in 1:nt) {
    # prior for exploitation rate
    U[t] ~ dbeta(1,1)
    
    # total catch and abundance
    C_tot[t] <- sum(C[t,1:ns])
    N_tot[t] <- sum(N[t,1:ns])
    
    # calculate age composition of aggregate run
    for (a in 1:na) {
      N_a_tot[t,a] <- sum(N_tas[t,a,1:ns])
      q_tot[t,a] <- N_a_tot[t,a]/N_tot[t]
    }
    
    # fit to aggregate harvest
    log_C_tot[t] <- log(C_tot[t])
    C_tot_t_obs[t] ~ dlnorm(log_C_tot[t],tau_C_obs[t])
  }
  
  # fit to population-specific escapement
  # use a vectorized format to handle any missing escapement data
  for (i in 1:S_obs_n) {
    log_S[i] <- log(S[S_obs_t[i],S_obs_s[i]])
    S_obs[i] ~ dlnorm(log_S[i],tau_S_obs[i])
  }
  
  # fit to age composition data: for years where data are not stock-specific (not_ss)
  for (t in 1:n_not_ss_years) {
    x_tas_obs[not_ss_years[t],1:na,1] ~ dmulti(q_tot[not_ss_years[t],1:na],ESS_ts[not_ss_years[t],1])
  }
  
  # fit to age composition data: for years where data are stock-specific (ss)
  for (t in 1:n_ss_years) {
    for (j in 1:ns) {
      x_tas_obs[ss_years[t],1:na,j+1] ~ dmulti(q_stock[ss_years[t],1:na,j],ESS_ts[ss_years[t],j+1])
    }
  }
}

# write the model to a text file
postpack::write_model(jags_model, jags_model_file)

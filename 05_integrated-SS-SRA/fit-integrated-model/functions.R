
# these functions were adapted from the code used to fit the integrated model to Kuskokwim River data
# much of that code was bundled into functions found in two packages:
# "FitSR" (https://www.github.com/bstaton1/FitSR)
# "SimSR" (https://www.github.com/bstaton1/SimSR)

##### RAW_DATA_PREP() #####
# ADAPTED FROM FitSR::kusko_data_prep()

# prepares raw data files into a form used by other functions

# S_dat = read.csv("inputs/esc-data.csv", stringsAsFactors = F)
# H_dat = read.csv("inputs/catch-data.csv", stringsAsFactors = F)
# age_dat = read.csv("inputs/age-data.csv", stringsAsFactors = F)

raw_data_prep = function (S_dat, H_dat, age_dat) {
  
  # prep escapement data
  S_tot_t_obs = tapply(S_dat$mean, S_dat$year, sum)
  S_dat$mean[S_dat$obs == 0] = NA
  years = unique(S_dat$year)
  nt = length(years)
  stocks = unique(S_dat$stock)
  ns = length(stocks)
  S_ts_obs = round(reshape2::dcast(S_dat, year ~ stock, value.var = "mean"))
  rownames(S_ts_obs) = years
  S_ts_obs = S_ts_obs[, -1]
  S_ts_obs = S_ts_obs[,stocks]
  cv_S_ts_obs = reshape2::dcast(S_dat, year ~ stock, value.var = "cv")
  rownames(cv_S_ts_obs) = years
  cv_S_ts_obs = cv_S_ts_obs[,-1]
  cv_S_ts_obs = cv_S_ts_obs[,stocks]
  sig_S_ts_obs = sqrt(log(cv_S_ts_obs^2 + 1))
  
  # prep catch data
  C_tot_t_obs = H_dat$mean
  cv_C_tot_t_obs = H_dat$cv
  sig_C_tot_t_obs = sqrt(log(cv_C_tot_t_obs^2 + 1))
  U_t_obs = C_tot_t_obs/(S_tot_t_obs + C_tot_t_obs)
  v = rep(1, ns)
  
  # dimensions
  a_min = 4
  a_max = 7
  na = a_max - a_min + 1
  ages = a_min:a_max
  ny = nt + na - 1
  age_stocks = which(stocks %in% age_dat$stock)
  n_age_stocks = length(age_stocks)
  
  # age data: ss means "stock-specific"
  ss_years = unique(age_dat[age_dat$stock != "aggregate" & !is.na(age_dat$n),"year"])
  not_ss_years = years[!(years %in% ss_years)]
  
  ss_year_ind = which(years %in% ss_years)
  not_ss_year_ind = which(years %in% not_ss_years)
  
  first_ss_age_yr = min(which(years %in% ss_years))
  last_ss_age_yr  = max(which(years %in% ss_years))
  
  x_tas_obs = array(NA, dim = c(nt, na, ns+1))
  
  jnames = unique(age_dat$stock)
  for (j in 1:(ns + 1)) {
    if (jnames[j] == "aggregate") {
      x_tas_obs[,,j] = as.matrix(age_dat[age_dat$stock == jnames[j],c("a4", "a5", "a6", "a7")])
    } else {
      x_tas_obs[ss_year_ind,,j] = as.matrix(age_dat[age_dat$stock == jnames[j] & age_dat$year %in% ss_years,c("a4", "a5", "a6", "a7")])
    }
  }
  
  # bundle
  obs = list(C_tot_t_obs = C_tot_t_obs, S_ts_obs = as.matrix(S_ts_obs), 
             x_tas_obs = x_tas_obs, sig_S_ts_obs = as.matrix(sig_S_ts_obs), 
             sig_C_t_obs = sig_C_tot_t_obs, U_t_obs = U_t_obs)
  
  params = list(ns = ns, nt = nt, a_min = a_min, a_max = a_max, 
                na = na, ages = a_min:a_max, ny = ny, stocks = stocks, 
                v = v, ss_years = ss_year_ind, not_ss_years = not_ss_year_ind)
  
  list(obs = obs, params = params)
}

##### gen_Rys_obs() #####
# ADAPTED FROM SimSR::gen_Rys_obs()

# reconstructs recruitment using standard brood table calculations

gen_Rys_obs = function (params, obs) {
  output = with(append(params, obs), {
    H_ts_obs = matrix(NA, nt, ns)
    R_ys_obs = matrix(NA, ny, ns)
    S_tas_obs = array(NA, dim = c(nt, na, ns))
    H_tas_obs = array(NA, dim = c(nt, na, ns))
    N_tas_obs = array(NA, dim = c(nt, na, ns))
    if (length(dim(x_tas_obs)) == 2) {
      q_tas_obs = t(apply(x_tas_obs, 1, function(x) x/sum(x)))
    } else {
      q_tas_obs = t(apply(x_tas_obs[,,1], 1, function(x) x/sum(x)))
    }
    
    for (s in 1:ns) {
      H_ts_obs[,s] = (S_ts_obs[,s] * U_t_obs)/(1 - U_t_obs)
      S_tas_obs[,,s] = apply(q_tas_obs, 2, function(x) x * S_ts_obs[,s])
      H_tas_obs[,,s] = apply(q_tas_obs, 2, function(x) x * H_ts_obs[,s])
      N_tas_obs[,,s] = S_tas_obs[,,s] + H_tas_obs[,,s]
      for (y in 1:ny) {
        if (y <= (nt - 3)) {
          brd.yr.runs = diag(N_tas_obs[y:(y + na - 1),,s])
          R_ys_obs[y + na - 1, s] = sum(brd.yr.runs, 
                                        na.rm = all(!is.na(brd.yr.runs)))
        }
        else {
          (next)()
        }
      }
    }
    N_ts_obs = apply(N_tas_obs, 3, rowSums)
    S_ind = 1:(nt - a_max)
    R_ind = (a_max + 1):(ny - na + 1)
    list(
      R_ys_obs = R_ys_obs, 
      n_S_obs = apply(S_ts_obs, 2, function(x) sum(!is.na(x))),
      n_R_obs = apply(R_ys_obs, 2, function(x) sum(!is.na(x))),
      n_SR_obs = sapply(1:ns, function(s) {sum(!is.na(S_ts_obs[S_ind, s]) & !is.na(R_ys_obs[R_ind,s]))}))
  })
  obs = append(obs, output)
  return(obs)
}

##### JAGS_DATA_PREP #####
# ADAPTED FROM FitSR::ssm_data_prep()

# prepares info to be passed to JAGS

jags_data_prep = function(params, obs) {
  output = with(append(params, obs), {
    S_ts_obs_m = S_ts_obs
    S_ts_obs_v = as.numeric(S_ts_obs_m)
    S_obs_s = rep(1:ns, each = nt)
    S_obs_t = rep(1:nt, ns)
    no_na_yrs = which(!is.na(S_ts_obs_v))
    S_obs = S_ts_obs_v[no_na_yrs]
    S_obs_s = S_obs_s[no_na_yrs]
    S_obs_t = S_obs_t[no_na_yrs]
    S_obs_n = length(S_obs)
    sig_S_obs_ts_v = as.numeric(sig_S_ts_obs)
    sig_S_obs = sig_S_obs_ts_v[no_na_yrs]
    x_tas_obs[is.na(x_tas_obs)] = 0
    ESS_ts = apply(x_tas_obs, 3, rowSums)
    
    out = list(ns = ns, nt = nt, ny = ny, na = na, a_max = a_max, 
         C_tot_t_obs = C_tot_t_obs, tau_C_obs = 1/sig_C_t_obs^2, 
         v = v, S_obs = S_obs, S_obs_t = S_obs_t, S_obs_s = S_obs_s, 
         S_obs_n = S_obs_n, tau_S_obs = 1/sig_S_obs^2, x_tas_obs = x_tas_obs, 
         ESS_ts = ESS_ts,
         R_wish = diag(rep(1, ns)), df_wish = ns + 1,
         ss_years = ss_years, not_ss_years = not_ss_years,
         n_ss_years = length(ss_years), n_not_ss_years = length(not_ss_years))
    out
  })
  return(output)
}

##### FUNCTIONS FOR GENERATING INITIAL VALUES FOR SSM ####

# ADAPTED FROM FitSR::lm_data_prep()

# prepares info for fitting a Ricker SR model using regression
# see Staton et al. (2020; 10.1139/cjfas-2019-0281) for details
# this model is fitted to get reasonable starting values for the MCMC algorithm

lm_data_prep = function (params, obs) {
  output = with(append(params, obs), {
    S_ind = 1:(nt - a_max)
    R_ind = (a_max + 1):(ny - na + 1)
    S_ys_reg = S_ts_obs[S_ind, ]
    R_ys_reg = R_ys_obs[R_ind, ]
    log_RPS_ys = log(R_ys_reg/S_ys_reg)
    lm_dat = data.frame(log_RPS = as.numeric(log_RPS_ys), 
                        S_ys = as.numeric(S_ys_reg), stock = rep(1:ns, each = nrow(S_ys_reg)))
    lm_dat = lm_dat[!is.na(lm_dat$log_RPS), ]
    list(ns = ns, nobs = nrow(lm_dat), obs_log_RPS_lme = lm_dat$log_RPS, 
         obs_log_RPS_lm = lm_dat$log_RPS, S_obs = lm_dat$S_ys, 
         stock = lm_dat$stock)
  })
  return(output)
}

# take estimates from regression model fit and get estimates of Umsy and Smsy
# there is no analytical solution to go from alpha to Umsy (unlike the reverse)
# so to obtain this quantity we need an approximation
# ADAPTED FROM FitSR::gen_lm_mgmt()

gen_lm_mgmt = function (alpha, beta) {
  log_alpha = log(alpha)
  U_msy = log_alpha * (0.5 - (0.65 * log_alpha^1.27)/(8.7 + log_alpha^1.27))
  U_msy[U_msy == "NaN"] = 0
  U_msy[U_msy < 0] = 0
  S_msy = U_msy/beta
  return(list(U_msy = U_msy, S_msy = S_msy))
}

# creates a list of initial values for the SSM
# ADAPTED FROM FitSR::gen_ssm_inits()
gen_ssm_inits = function (params, obs, n_chains) {
  output = with(append(params, obs), {
    lm_data = lm_data_prep(params = params, obs = obs)
    lm_alpha = NULL
    lm_beta = NULL
    for (s in 1:ns) {
      tmp_S = lm_data$S_obs[lm_data$stock == s]
      tmp_log_RPS = lm_data$obs_log_RPS_lm[lm_data$stock == s]
      tmp_fit = tryCatch({
        lm(tmp_log_RPS ~ tmp_S)
      }, error = function(e) NULL)
      if (is.null(tmp_fit)) {
        lm_alpha = c(lm_alpha, NA)
        lm_beta = c(lm_beta, NA)
      }
      else {
        lm_alpha = c(lm_alpha, unname(exp(coef(tmp_fit)[1])))
        lm_beta = c(lm_beta, unname(abs(coef(tmp_fit)[2])))
      }
    }
    lm_alpha[lm_alpha <= 1] = 1.5
    lm_alpha[lm_alpha > 20] = 20
    lm_alpha[is.na(lm_alpha)] = mean(lm_alpha, na.rm = T)
    lm_beta[is.na(lm_beta)] = mean(lm_beta, na.rm = T)
    lm_mgmt = gen_lm_mgmt(lm_alpha, lm_beta)
    inits = list()
    for (i in 1:n_chains) {
      inits[[i]] = list(
        U_msy = {
          y = rbeta(ns, lm_mgmt$U_msy * 100, (1 - lm_mgmt$U_msy) * 100)
          ifelse(y < 0.1, 0.15, y)
        }, 
        log_S_msy = log(rlnorm(ns, log(lm_mgmt$S_msy), 0.1)), 
        log_R = apply(R_ys_obs, 2, function(x) {
          mu = mean(x, na.rm = T)
          x[is.na(x)] = mu
          log(rlnorm(length(x), log(x), 0.2))
        }),
        D_scale = runif(1, 0.08, 0.12),
        phi = runif(params$ns, 0, 0.5))
    }
    inits
  })
  return(output)
}



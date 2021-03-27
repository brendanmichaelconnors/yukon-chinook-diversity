
##### SESSION INITIALIZATION #####

# SET THE WORKING DIRECTORY TO THE LOCATION OF THIS FILE

# clear the workspace
rm(list = ls(all = T))

# MCMC settings
short_run = F  # ~1.5 minutes -- for testing out this code. Totally inadequate for inference
med_run = T    # ~1.0 hours -- for testing plotting code. Quasi-valid for inference
long_run = F   # ~1.5 days -- used in manuscript analysis. Passes all diagnostic checks

# read in other misc functions for this analysis
source("functions.R")

# export the JAGS model to a text file
source("ssm-model.R")

# path to inputs and file names
in_dir = "inputs"
C_file = file.path(in_dir, "catch-data.csv")
S_file = file.path(in_dir, "esc-data.csv")
A_file = file.path(in_dir, "age-data.csv")

# path to outputs and file names
out_dir = "../posterior-samples"; if (!dir.exists(out_dir)) dir.create(out_dir)
post_file = "integrated-posterior.rds"

##### DATA PREPARATIONS #####
inputs = raw_data_prep(
  S_dat = read.csv(S_file, stringsAsFactors = F),
  H_dat = read.csv(C_file, stringsAsFactors = F),
  age_dat = read.csv(A_file, stringsAsFactors = F)
)

# extract the observed data
obs = inputs$obs

# extract various dimensional quantities
params = inputs$params

# append brood-table reconstructed recruitments to the observations
# useful for plotting and generating initial values
obs = gen_Rys_obs(params, obs)

# create the data list to pass to JAGS
jags_data = jags_data_prep(params, obs)

##### PARAMETERS TO MONITOR #####
jags_params = c(
  "alpha", "beta", "U_msy", "S_msy", "U", "pi", "R_eq",
  "sigma_R", "rho_mat", "phi",  "S", "p",
  "R", "log_resid", "Sigma_R", "q", "q_tot", "q_stock", "C_tot", "D_sum"
)

##### MCMC DIMENSIONS #####
if (sum(c(short_run, med_run, long_run) != 1)) {
  stop ("only one of 'short_run', 'med_run', or 'long_run' must be TRUE")
}

if (short_run) jags_dims = c(ni = 500, nb = 1, nt = 1, nc = 4, na = 100)
if (med_run) jags_dims = c(ni = 20000, nb = 15000, nt = 10, nc = 4, na = 1000)
if (long_run) jags_dims = c(ni = 5e5, nb = 1e6, nt = 200, nc = 4, na = 50000)

##### MCMC INITIAL VALUES #####
jags_inits = gen_ssm_inits(params, obs, jags_dims["nc"])

##### RUN THE SAMPLER #####
start = Sys.time(); cat("MCMC Started At:", format(start), "\n")
# run the sampler
post = jagsUI::jags.basic(
  data = jags_data,
  inits = jags_inits,
  parameters.to.save = jags_params,
  model.file = jags_model_file,
  n.chains = jags_dims["nc"],
  n.adapt = jags_dims["na"],
  n.iter = sum(jags_dims[c("ni", "nb")]),
  n.burnin = jags_dims["nb"],
  n.thin = jags_dims["nt"],
  parallel = T,
  verbose = T,
  factories = "bugs::MNormal sampler FALSE",
  save.model = F
)
stop = Sys.time()

# delete the JAGS model file
unlink(jags_model_file)

##### SAVE THE OUTPUT #####
saveRDS(post, file = file.path(out_dir, post_file))
stop - start

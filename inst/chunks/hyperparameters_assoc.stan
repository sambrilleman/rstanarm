  // hyperparameter values are set to 0 if there is no prior
  vector[e_A]          e_prior_mean_for_assoc;
  vector<lower=0>[e_A] e_prior_scale_for_assoc;
  vector<lower=0>[e_A] e_prior_df_for_assoc;
  real<lower=0>        e_global_prior_scale_for_assoc; // for hs priors only
  real<lower=0>        e_global_prior_df_for_assoc;  

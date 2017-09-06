  vector[e_A] e_z_alpha; // primitive assoc params

  // parameters for priors on assoc params 
  real<lower=0> e_global_for_assoc[e_hs_for_assoc];
  vector<lower=0>[(e_hs_for_assoc>0)*e_A] e_local_for_assoc[e_hs_for_assoc];
  vector<lower=0>[e_A] e_mix_for_assoc[e_prior_dist_for_assoc == 5 || e_prior_dist_for_assoc == 6];
  real<lower=0> e_ool_for_assoc[e_prior_dist_for_assoc == 6];

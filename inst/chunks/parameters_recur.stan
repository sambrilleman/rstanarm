int<lower=0,upper=1> has_recurrent; // 1 = has recurrent event submodel
real r_gamma[r_has_intercept]; // intercept for event submodel 
vector[r_K] r_z_beta;          // primitive log hazard ratios
real e_z_fbeta[has_recurrent]; // coef to scale frailty term in event submodel 

// unscaled basehaz params, either:
//   - weibull shape parameter
//   - b-spline coefs on log basehaz
//   - coefs for piecewise constant basehaz
vector<lower=(basehaz_type == 1 ? 0 : negative_infinity())>[basehaz_df*has_recurrent] r_aux_unscaled;       

// parameters for priors on log haz ratios
real<lower=0> r_global[r_hs*has_recurrent];
vector<lower=0>[(r_hs>0)*r_K] r_local[r_hs*has_recurrent];
vector<lower=0>[r_K] r_mix[(r_prior_dist == 5 || r_prior_dist == 6)*has_recurrent];
real<lower=0> r_ool[(r_prior_dist == 6)*has_recurrent];

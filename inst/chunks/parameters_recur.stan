real r_gamma[r_has_intercept]; // intercept for event submodel 
vector[r_K] r_z_beta;          // primitive log hazard ratios
vector[has_recurrent] e_z_lambda; // coef to scale frailty term in event submodel 
vector[Npat*has_recurrent] z_frailty; // standardised frailty terms
vector<lower=0>[has_recurrent] frailty_aux;

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

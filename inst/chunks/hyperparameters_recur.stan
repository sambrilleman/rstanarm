// hyperparameter values are set to 0 if there is no prior
vector[1]                   e_prior_mean_for_lambda;
vector[r_K]                 r_prior_mean;
real                        r_prior_mean_for_intercept;
vector[basehaz_df]          r_prior_mean_for_aux;
vector<lower=0>[1]          e_prior_scale_for_lambda;
vector<lower=0>[r_K]        r_prior_scale;
real<lower=0>               r_prior_scale_for_intercept;
vector<lower=0>[basehaz_df] r_prior_scale_for_aux;
vector<lower=0>[1]          e_prior_df_for_lambda;
vector<lower=0>[r_K]        r_prior_df;
real<lower=0>               r_prior_df_for_intercept; 
vector<lower=0>[basehaz_df] r_prior_df_for_aux; 
real<lower=0>               r_global_prior_scale; // for hs priors only
real<lower=0>               r_global_prior_df;  
real<lower=0>               f_prior_scale_for_aux;

// prior family: 0 = none, 1 = log-normal, 2 = gamma
int<lower=0,upper=2> frailty_dist;

// prior family: 0 = none, 1 = normal, 2 = student_t, 3 = hs, 4 = hs_plus, 
//   5 = laplace, 6 = lasso
int<lower=0,upper=3> e_prior_dist_for_lambda;
int<lower=0,upper=6> r_prior_dist;
int<lower=0,upper=2> r_prior_dist_for_intercept;

// prior family: 0 = none, 1 = normal, 2 = student_t, 3 = exponential
int<lower=0,upper=3> r_prior_dist_for_aux; // prior for basehaz params

// data for recurrent event submodel
int<lower=0,upper=1> has_recurrent; // 1 = has recurrent event submodel
int<lower=0> r_K;           // num. of predictors in event submodel
int<lower=0,upper=1> r_has_intercept; // 1 = yes
int<lower=0> nrow_r_Xq;     // num. rows in event predictor matrix at quad points
matrix[nrow_r_Xq,r_K] r_Xq; // predictor matrix (event submodel) at quadpoints, centred
vector[nrow_r_Xq] r_times;  // event times and unstandardised quadrature points
matrix[nrow_r_Xq,basehaz_df*has_recurrent] r_basehaz_X; // design matrix (basis terms) for baseline hazard
vector[nrow_r_Xq] r_d;      // event indicator, followed by dummy indicator for quadpoints
vector[r_K] r_xbar;         // predictor means (event submodel)
int<lower=0> Nri[Npat*has_recurrent];     // num. of recurrent events for each individual
int<lower=0> Nrtotal;       // total num. of recurrent events
int<lower=0,upper=Npat> e_uindices[nrow_e_Xq*has_recurrent];  // indexing for frailty terms in recurrent event submodel
int<lower=0,upper=Npat> r_uindices[nrow_r_Xq*has_recurrent];  // indexing for frailty terms in recurrent event submodel

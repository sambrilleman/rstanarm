#include "Columbia_copyright.stan"
#include "Brilleman_copyright.stan"
#include "license.stan" // GPL3+

// Shared parameter joint model
functions {
  #include "common_functions.stan"
  #include "bernoulli_likelihoods.stan"
  #include "binomial_likelihoods.stan"
  #include "continuous_likelihoods.stan"
  #include "count_likelihoods.stan"
  #include "jm_functions.stan"
}
data {
  #include "NKX.stan" // declares N, K, X, xbar, dense_X, nnz_x, w_x, v_x, u_x  
  
  // declares M, N{M,_real,_int,01}, idx{_real,_int,_K,_hs2,hs4}, KM,
  //   sum_has_intercept{_nob,_lob,_upb}, sum_has_aux
  #include "dimensions_mvglm.stan" 
  
  // declares prior_PD, has_intercept{_nob,_lob,_upb}, family, link, prior_dist, 
  //   prior_dist_for_intercept, prior_dist_for_aux
  #include "data_mvglm.stan"     // same as data_glm.stan, but arrays of size M 
  #include "data2_mvglm.stan"    // declares y_{real,int}, has_aux 
  #include "weights_offset.stan" // declares has_weights, weights, has_offset, offset
 
  // declares t, p[t], l[t], q, len_theta_L, shape, scale, 
  //   {len_}concentration, {len_}regularization
  #include "glmer_stuff.stan"  
  #include "glmer_stuff2.stan" // declares num_not_zero, w, v, u, special_case
  #include "mvmer_stuff.stan"  // declares pmat, qmat, q1, q2
  int<lower=1,upper=t> t_i;    // index of grouping factor corresponding to patient-level

  // declares e_prior_dist{_for_intercept,_for_aux}, Npat{_times_}quadnodes, quadweight, 
  //   basehaz_{type,df,X}, nrow_e_Xq, e_{K,Xq,times,d,xbar,weights,weights_rep}  
  #include "data_event.stan"

  // declares e_prior_dist_for_{frailty,frscale}, r_prior_dist{_for_intercept,_for_aux}, 
  //   r_basehaz_X, nrow_r_Xq, r_{K,Xq,times,d,xbar}  
  #include "data_recur.stan"

  // declares e_A, e_prior_dist_for_assoc, assoc, assoc_uses, has_assoc, {sum_}size_which_b, 
  //   which_b_zindex, {sum_}size_which_coef, which_coef_{zindex,xindex}, 
  //   {sum_}a_K_data, {sum_,sum_size_}which_interactions, y_Xq_{eta,eps,lag,auc,data},
  //   {nnz,w,v,u}_Zq_{eta,eps,lag,auc}, nrow_y_Xq, nrow_y_Xq_auc, 
  //   auc_quadnodes, auc_quadweights   
  #include "data_assoc.stan"
  
  // declares {e_,r_}prior_{mean, scale, df}, {e_,r_}prior_{mean, scale, df}_for_intercept, 
  //   {e_,r_}prior_{mean, scale, df}_for_aux, e_prior_{mean, scale, df}_for_{assoc,frailty,frscale},
  //   {e_,r_}global_prior_{df,scale}
  #include "hyperparameters_mvglm.stan" // same as hyperparameters.stan, but arrays of size M
  #include "hyperparameters_event.stan" 
  #include "hyperparameters_assoc.stan" 
  #include "hyperparameters_recur.stan" 
}
transformed data {
  int<lower=0> e_hs = get_nvars_for_hs(e_prior_dist);                 
  int<lower=0> e_hs_for_assoc = get_nvars_for_hs(e_prior_dist_for_assoc);                 
  int<lower=0> r_hs = get_nvars_for_hs(r_prior_dist);                 
  int<lower=1> V[special_case ? t : 0, N] = make_V(N, special_case ? t : 0, v);  // not used
  
  // declares poisson_max, hsM, idx_{global,local2,local4,mix,ool,noise}, 
  //   len_{global,local2,local4,mix,ool,noise}, {sqrt_,log_,sum_log_}y, 
  //   len_z_T, len_var_group, delta, is_continuous, pos, beta_smooth
  #include "tdata_mvglm.stan" 

  // defines hsM, idx_{global,local2,local4,mix,ool,noise}, 
  //   len_{global,local2,local4,mix,ool}, {sqrt_,log_,sum_log_}y, 
  //   len_z_T, len_var_group, delta, is_continuous, pos  
  #include "tdata2_mvglm.stan" 
}
parameters {
  // declares gamma_{nob,lob,upb}, z_beta, global, local{2,4}, mix, 
  //   ool, noise, aux_unscaled, z_b, z_T, rho, zeta, tau
  #include "parameters_mvglm.stan"
  // declares e_{gamma,z_beta,aux_unscaled,global,local,mix,ool}  
  #include "parameters_event.stan"
  // declares e_z_alpha, e_{global,local,mix,ool}_for_assoc
  #include "parameters_assoc.stan"  
  // declares r_{z_beta,global,local,mix,ool}, e_z_fbeta
  #include "parameters_recur.stan" 
}
transformed parameters { 
  // parameters for event submodel
  vector[e_K] e_beta;               // log hazard ratios
  vector[e_A] e_alpha;              // assoc params
  vector[basehaz_df] e_aux;         // basehaz params  
  vector[has_recurrent] e_fbeta;    // coef on frailty term in event submodel
  vector[r_K] r_beta;               // log hazard ratios
  vector[basehaz_df*has_recurrent] r_aux; // basehaz params 
  vector[Npat*has_recurrent] frailty; // frailty terms
  #include "tparameters_mvglm.stan" // defines aux, beta, b{_not_by_model}, theta_L
  e_beta = generate_beta(e_z_beta, e_prior_dist, e_prior_mean, 
                         e_prior_scale, e_prior_df, e_global, e_local,
                         e_global_prior_scale, e_ool, e_mix);  
  e_alpha = generate_beta(e_z_alpha, e_prior_dist_for_assoc, e_prior_mean_for_assoc, 
                          e_prior_scale_for_assoc, e_prior_df_for_assoc, 
                          e_global_for_assoc, e_local_for_assoc,
                          e_global_prior_scale_for_assoc, e_ool_for_assoc, 
                          e_mix_for_assoc);         
  e_aux  = generate_aux(e_aux_unscaled, e_prior_dist_for_aux,
                        e_prior_mean_for_aux, e_prior_scale_for_aux);
  if (has_recurrent == 1) {
    r_beta = generate_beta(r_z_beta, r_prior_dist, r_prior_mean, 
                           r_prior_scale, r_prior_df, r_global, r_local,
                           r_global_prior_scale, r_ool, r_mix);  
    r_aux  = generate_aux(r_aux_unscaled, r_prior_dist_for_aux,
                          r_prior_mean_for_aux, r_prior_scale_for_aux);
    e_fbeta = generate_aux(e_z_fbeta, e_prior_dist_for_frscale, e_prior_mean_for_frscale, 
                           e_prior_scale_for_frscale);  
    frailty = generate_aux(z_frailty, 1, frailty_mean, frailty_sd);  
  }
  if (t > 0) {
    theta_L = make_theta_L(len_theta_L, p, 1.0, tau, scale, zeta, rho, z_T);
    b_not_by_model = make_b(z_b, theta_L, p, l);
    if (M == 1) b = b_not_by_model;
	  else b = reorder_b(b_not_by_model, p, pmat, q1, q2, qmat, l, M);
  }
}
model {
  vector[nrow_e_Xq] e_eta_q; // eta for event submodel (at event and quad times)  
  vector[nrow_r_Xq*has_recurrent] r_eta_q; // eta for recurrent event submodel (at event and quad times)  

  //---- Log-lik for longitudinal submodels
  
  int aux_mark = 1; 
  #include "make_eta.stan" // defines eta
  if (t > 0) {
    #include "eta_add_Zb.stan"
  }
  for (m in 1:M) {
    vector[NM[m]] eta_tmp;	            // eta for just one submodel 
    eta_tmp = eta[idx[m,1]:idx[m,2]];   // eta for just one submodel 
    #include "eta_intercept_mvmer.stan"	// adds intercept or shifts eta
    if (family[m] == 8) {  // poisson-gamma mixture
	  #include "eta_add_noise_mvmer.stan"
    }    
    #include "mvmer_lp.stan" // increments target with long log-liks
    if (has_aux[m] == 1) aux_mark = aux_mark + 1;
  }
  
  //----- Log-lik for event submodel (GK quadrature)
  
  // Event submodel: linear predictor at event and quad times
  if (e_K > 0) e_eta_q = e_Xq * e_beta;
  else e_eta_q = rep_vector(0.0, nrow_e_Xq);
  if (e_has_intercept == 1) e_eta_q = e_eta_q + e_gamma[1];
  else e_eta_q = e_eta_q + dot_product(e_xbar, e_beta);     
  if (has_recurrent == 1) {
    e_eta_q = e_eta_q + e_fbeta[1] * frailty[e_uindices];
  }
  if (assoc == 1) { 
    // declares y_eta_q{_eps,_lag,_auc}, y_eta_qwide{_eps,_lag,_auc}, 
	  //   y_q_wide{_eps,_lag,_auc}, mark{2,3}
    #include "assoc_definitions.stan"  
    #include "assoc_prepwork.stan"
    #include "assoc_evaluate.stan"
  }
  { 
    // declares log_basehaz, ll_{haz_q,haz_eventtime,surv_eventtime,event}
	  #include "event_lp.stan" // increments target with event log-lik
  }
  
  //----- Log-lik for recurrent event submodel (GK quadrature)
  
  if (has_recurrent == 1) {
    // Recurrent event submodel: linear predictor at event and quad times
    if (r_K > 0) r_eta_q = r_Xq * r_beta;
    else r_eta_q = rep_vector(0.0, nrow_r_Xq);
    if (r_has_intercept == 1) r_eta_q = r_eta_q + r_gamma[1];
    else r_eta_q = r_eta_q + dot_product(r_xbar, r_beta); 
    r_eta_q = r_eta_q + frailty[r_uindices];
    { 
      // declares log_basehaz, ll_{haz_q,haz_eventtime,surv_eventtime,event}
  	  #include "recur_lp.stan" // increments target with event log-lik
    }  
  }
  
  //----- Log-priors

  // increment target with priors for longitudinal submodel params
  for (m in 1:M) { 
    if (has_aux[m] == 1) {
      aux_mark = sum(has_aux[1:m]);
      aux_lp(aux_unscaled[aux_mark], prior_dist_for_aux[m], prior_scale_for_aux[m], prior_df_for_aux[m])
    } 
  }
  #include "priors_mvglm.stan"  
  
  // increment target with priors for event submodel params
  beta_lp(e_z_beta, e_prior_dist, e_prior_scale, e_prior_df, 
          e_global_prior_df, e_local, e_global, e_mix, e_ool)
  beta_lp(e_z_alpha, e_prior_dist_for_assoc, e_prior_scale_for_assoc, 
          e_prior_df_for_assoc, e_global_prior_df_for_assoc, e_local_for_assoc, 
          e_global_for_assoc, e_mix_for_assoc, e_ool_for_assoc)
  basehaz_lp(e_aux_unscaled, e_prior_dist_for_aux, 
             e_prior_scale_for_aux, e_prior_df_for_aux)
  if (e_has_intercept == 1) 
    gamma_lp(e_gamma[1], e_prior_dist_for_intercept, e_prior_mean_for_intercept, 
             e_prior_scale_for_intercept, e_prior_df_for_intercept);  
  
  // increment target with priors for recurrent event submodel params
  if (has_recurrent == 1) {
    gamma_lp(e_z_fbeta[1], e_prior_dist_for_frscale, e_prior_mean_for_frscale[1],
             e_prior_scale_for_frscale[1], e_prior_df_for_frscale[1])
    beta_lp(r_z_beta, r_prior_dist, r_prior_scale, r_prior_df, 
            r_global_prior_df, r_local, r_global, r_mix, r_ool)
    basehaz_lp(r_aux_unscaled, r_prior_dist_for_aux,
               r_prior_scale_for_aux, r_prior_df_for_aux)
    if (r_has_intercept == 1) 
      gamma_lp(r_gamma[1], r_prior_dist_for_intercept, r_prior_mean_for_intercept, 
               r_prior_scale_for_intercept, r_prior_df_for_intercept);
    // frailty terms
    target += normal_lpdf(z_frailty | 0, 1);
  }

  // increment target with priors for group-specific params
  if (t > 0) decov_lp(z_b, z_T, rho, zeta, tau, 
                      regularization, delta, shape, t, p);
}
generated quantities {
  #include "gen_quantities_mvmer.stan"  // defines alpha, mean_PPD
}

vector[nrow_r_Xq] log_basehaz;  // baseline hazard evaluated at quadpoints
vector[nrow_r_Xq] ll_haz_q;     // log hazard contribution to the log likelihood for the event model at recurrent event time and quad points
vector[Nrtotal] ll_haz_recurtime;  // log hazard contribution to the log likelihood for the event model AT the recurrent event time only
vector[Npat_times_quadnodes] ll_surv_eventtime; // log survival contribution to the log likelihood for the event model AT the event time
real ll_recur;                  // log likelihood for the recurrent event model  

// Log baseline hazard at recurrent event and quad times
if (basehaz_type == 1) 
  log_basehaz = log(r_aux[1]) + r_basehaz_X * (r_aux - 1);
else log_basehaz = r_basehaz_X * r_aux;	

// Log hazard at event and quad times
ll_haz_q = r_d .* (log_basehaz + r_eta_q);

// Log hazard contribution to the likelihood
ll_haz_recurtime = segment(ll_haz_q, 1, Nrtotal);

// Log survival contribution to the likelihood (by summing over the 
// quadrature points to get the approximate integral)
// NB quadweight already incorporates the (b-a)/2 scaling such that the
// integral is evaluated over limits (a,b) rather than (-1,+1)
ll_surv_eventtime = quadweight .* 
  exp(segment(ll_haz_q, (Nrtotal + 1), Npat_times_quadnodes));        

// Log likelihood for event model
if (has_weights == 0) { // unweighted log likelihood
  ll_recur = sum(ll_haz_recurtime) - sum(ll_surv_eventtime);
} 
else {  // weighted log likelihood
  // not yet implemented
}				    

// Log-likelihoods for event submodel  
if (prior_PD == 0) target += ll_recur;  

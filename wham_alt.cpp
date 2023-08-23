#include <TMB.hpp>
#include <iostream>

#include "helper_functions.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(n_years_catch);
  DATA_INTEGER(n_years_indices);
  DATA_INTEGER(n_years_model);
  DATA_INTEGER(n_ages);
  DATA_INTEGER(n_fleets);
  DATA_INTEGER(n_indices);
  DATA_INTEGER(n_selblocks);
  DATA_IVECTOR(selblock_models);
  DATA_IMATRIX(selblock_pointer_fleets);
  DATA_IMATRIX(selblock_pointer_indices);
  DATA_IVECTOR(age_comp_model_fleets);
  DATA_IVECTOR(n_age_comp_pars_fleets);
  DATA_IVECTOR(age_comp_model_indices);
  DATA_IVECTOR(n_age_comp_pars_indices);
  DATA_VECTOR(fracyr_SSB);
  DATA_MATRIX(mature);
  DATA_IVECTOR(waa_pointer_fleets);
  DATA_INTEGER(waa_pointer_totcatch);
  DATA_IVECTOR(waa_pointer_indices);
  DATA_INTEGER(waa_pointer_ssb);
  DATA_INTEGER(waa_pointer_jan1);
  DATA_ARRAY(waa);
  DATA_MATRIX(agg_catch);
  DATA_MATRIX(agg_catch_sigma);
  DATA_ARRAY(catch_paa); //n_fleets x n_years x n_ages
  DATA_IMATRIX(use_catch_paa);
  DATA_MATRIX(catch_Neff);
  DATA_IMATRIX(catch_aref);
  DATA_IVECTOR(units_indices);
  DATA_MATRIX(fracyr_indices);
  DATA_MATRIX(agg_indices);
  DATA_IMATRIX(use_indices);
  DATA_MATRIX(agg_index_sigma);
  DATA_IVECTOR(units_index_paa);
  DATA_ARRAY(index_paa); //n_indices x n_years x n_ages
  DATA_IMATRIX(use_index_paa);
  DATA_MATRIX(index_Neff);
  DATA_IMATRIX(index_aref);
  DATA_VECTOR(q_lower);
  DATA_VECTOR(q_upper);
  DATA_MATRIX(selpars_lower);
  DATA_MATRIX(selpars_upper);
  DATA_INTEGER(n_NAA_sigma);
  DATA_IVECTOR(NAA_sigma_pointers);
  DATA_INTEGER(recruit_model);
  DATA_INTEGER(n_M_re);
  DATA_IVECTOR(MAA_pointer); //n_ages
  DATA_IVECTOR(M_sigma_par_pointer); //n_M_re
  DATA_INTEGER(M_model); //0: just age-specific M, 1: Lorenzen M decline with age/weight
  DATA_INTEGER(N1_model); //0: just age-specific numbers at age, 1: 2 pars: log_N_{1,1}, log_F0, age-structure defined by equilibrium NAA calculations
  DATA_INTEGER(use_M_re);
  DATA_INTEGER(use_NAA_re);
  DATA_INTEGER(use_b_prior);
  DATA_INTEGER(random_recruitment);
  DATA_INTEGER(which_F_age); //which age of F to use for full total F for msy/ypr calculations
  DATA_INTEGER(use_steepness); // which parameterization to use for BH/Ricker S-R, if needed.
  DATA_INTEGER(bias_correct_pe); //bias correct lognormal process error?
  DATA_INTEGER(bias_correct_oe); //bias correct lognormal observation error?
  DATA_IVECTOR(Fbar_ages);
  DATA_INTEGER(simulate_state); //if 1 then state parameters will be simulated
  DATA_SCALAR(percentSPR); //percentage to use for SPR-based reference points
  DATA_INTEGER(steepness_year); //year of SPR0 to use for S-R fit
  DATA_INTEGER(R0_is_S0); //whether to use R0 or S0 in steepness formulation
  DATA_SCALAR(FMSY_startval) //which F to start at for newton method to get FMSY
  
  PARAMETER_VECTOR(mean_rec_pars);
  PARAMETER_VECTOR(logit_q);
  PARAMETER_VECTOR(log_F1);
  PARAMETER_MATRIX(F_devs);
  PARAMETER_VECTOR(log_N1_pars); //length = n_ages or 2
  PARAMETER_VECTOR(log_NAA_sigma);
  PARAMETER_MATRIX(logit_selpars); //n_selblocks x n_ages + 6 (n_ages for by age, 2 for logistic, 4 for double-logistic)
  PARAMETER_VECTOR(catch_paa_pars);
  PARAMETER_VECTOR(index_paa_pars);
  PARAMETER_MATRIX(log_NAA);
  PARAMETER_VECTOR(M_pars1);
  PARAMETER_VECTOR(M_sigma_pars); //up to n_M_re
  PARAMETER_MATRIX(M_re); //n_years-1 x n_M_re
  PARAMETER(log_b);
  PARAMETER_VECTOR(log_R); //n_years-1, if used.
  PARAMETER(log_R_sigma);
  PARAMETER_VECTOR(log_catch_sig_scale) //n_fleets
  PARAMETER_VECTOR(log_index_sig_scale) //n_indices

  Type nll= 0.0; //negative log-likelihood
  vector<int> any_index_age_comp(n_indices);
  vector<int> any_fleet_age_comp(n_fleets);
  vector<Type> SSB(n_years_model);
  matrix<Type> F(n_years_model,n_fleets);
  matrix<Type> log_F(n_years_model,n_fleets);
  array<Type> pred_CAA(n_years_model,n_fleets,n_ages);
  array<Type> pred_catch_paa(n_years_model,n_fleets,n_ages);
  matrix<Type> pred_catch(n_years_model,n_fleets);
  matrix<Type> log_pred_catch(n_years_model,n_fleets);
  array<Type> pred_IAA(n_years_model,n_indices,n_ages);
  array<Type> pred_index_paa(n_years_model,n_indices,n_ages);
  matrix<Type> pred_indices(n_years_model,n_indices);
  matrix<Type> NAA(n_years_model,n_ages);
  matrix<Type> pred_NAA(n_years_model,n_ages);
  array<Type> FAA(n_years_model,n_fleets,n_ages);
  array<Type> log_FAA(n_years_model,n_fleets,n_ages);
  matrix<Type> FAA_tot(n_years_model,n_ages);
  matrix<Type> ZAA(n_years_model,n_ages);
  array<Type> QAA(n_years_model,n_indices,n_ages);
  matrix<Type> selblocks(n_selblocks,n_ages);
  vector<Type> q(n_indices);
  vector<Type> t_paa(n_ages), t_pred_paa(n_ages);
  matrix<Type> selpars(n_selblocks,n_ages+6);
  //Type SR_a, SR_b, SR_R0, SR_h;
  for(int i = 0; i < n_selblocks; i++) for(int j = 0; j < n_ages + 6; j++) 
  {
    selpars(i,j) = selpars_lower(i,j) +  (selpars_upper(i,j)-selpars_lower(i,j))/(1.0+exp(-logit_selpars(i,j)));
  }
  REPORT(selpars);
  selblocks = get_selblocks(n_ages, n_selblocks, selpars, selblock_models);

  for(int i = 0; i < n_indices; i++)
  {
    any_index_age_comp(i) = 0;
    for(int y = 0; y < n_years_indices; y++) if(use_index_paa(y,i) == 1) any_index_age_comp(i) = 1;
  }
  for(int i = 0; i < n_fleets; i++)
  {
    any_fleet_age_comp(i) = 0;
    for(int y = 0; y < n_years_catch; y++) if(use_catch_paa(y,i) == 1) any_fleet_age_comp(i) = 1;
  }
  vector<Type> sigma2_log_NAA = exp(log_NAA_sigma*2.0);
 
  if(use_M_re == 1) 
  {
    vector<Type> sigma_M(n_M_re);
    for(int i = 0; i < n_M_re; i++) sigma_M(i) = exp(M_sigma_pars(M_sigma_par_pointer(i)-1));
    matrix<Type> nll_M(n_years_model-1, n_M_re);
    nll_M.setZero();
    
    if(M_model == 0) for(int i = 0; i < n_M_re; i++) 
    {
      Type mu = M_pars1(i);
      if(bias_correct_pe == 1) mu -= 0.5*exp(2*log(sigma_M(i)));
      nll_M(0,i) -= dnorm(M_re(0,i), mu, sigma_M(i), 1);
      SIMULATE if(simulate_state == 1) M_re(0,i) = rnorm(mu, sigma_M(i));
      for(int y = 1; y < n_years_model - 1; y++)  
      {
        mu = M_re(y-1,i);
        if(bias_correct_pe == 1) mu -= 0.5*exp(2*log(sigma_M(i)));
        nll_M(y,i) -= dnorm(M_re(y,i), mu, sigma_M(i), 1);
        SIMULATE if(simulate_state == 1) M_re(y,i) = rnorm(mu, sigma_M(i));
      }
    }
    else 
    {
      for(int i = 0; i < n_M_re; i++) 
      {
        Type mu = M_pars1(i);
        if(bias_correct_pe == 1) mu -= 0.5*exp(2*log(sigma_M(i)));
        nll_M(0,i) -= dnorm(M_re(0,i), mu, sigma_M(i), 1);
        SIMULATE if(simulate_state == 1) M_re(0,i) = rnorm(mu, sigma_M(i));
      }
      for(int i = 0; i < n_M_re; i++) for(int y = 1; y < n_years_model - 1; y++)  
      {
        Type mu = M_re(y-1,i);
        if(bias_correct_pe == 1) mu -= 0.5*exp(2*log(sigma_M(i)));
        nll_M(y,i) -= dnorm(M_re(y,i), mu, sigma_M(i), 1);
        SIMULATE if(simulate_state == 1) M_re(y,i) = rnorm(mu, sigma_M(i));
      }
    }
    SIMULATE REPORT(M_re);
    REPORT(nll_M);
    nll += nll_M.sum();
  }
  //see(nll);

  matrix<Type> MAA(n_years_model,n_ages);
  for(int i = 0; i < n_ages; i++) 
  {
    if(M_model == 0) //M by age
    {
      MAA(0,i) = exp(M_pars1(MAA_pointer(i)-1));
      for(int y = 1; y < n_years_model; y++) MAA(y,i) = exp(M_re(y-1,MAA_pointer(i)-1));
    }
    else //M_model == 1, allometric function of weight
    {
      MAA(0,i) = exp(M_pars1(0) - exp(log_b) * log(waa(waa_pointer_jan1-1,0,i)));
      if(use_M_re == 1) for(int y = 1; y < n_years_model; y++) MAA(y,i) = exp(M_re(y-1,0) - exp(log_b) * log(waa(waa_pointer_jan1-1,y,i)));
      else for(int y = 1; y < n_years_model; y++) MAA(y,i) = exp(M_pars1(0) - exp(log_b) * log(waa(waa_pointer_jan1-1,y,i)));
    }
  }

  for(int i = 0; i < n_indices; i++)
  {
    q(i) = q_lower(i) + (q_upper(i) - q_lower(i))/(1 + exp(-logit_q(i)));
    for(int y = 0; y < n_years_model; y++) 
    {
      for(int a = 0; a < n_ages; a++) QAA(y,i,a) = q(i) * selblocks(selblock_pointer_indices(y,i)-1,a);
    }
  }
  FAA_tot.setZero();
  for(int f = 0; f < n_fleets; f++)
  {
    log_F(0,f) = log_F1(f);
    F(0,f) = exp(log_F(0,f));
    for(int a = 0; a < n_ages; a++) 
    {
      FAA(0,f,a) = F(0,f) * selblocks(selblock_pointer_fleets(0,f)-1,a);
      log_FAA(0,f,a) = log(FAA(0,f,a));
      FAA_tot(0,a) = FAA_tot(0,a) + FAA(0,f,a);
    }
    for(int y = 1; y < n_years_model; y++) 
    {
      log_F(y,f) = log_F(y-1,f) + F_devs(y-1,f);
      F(y,f) = exp(log_F(y,f));
      for(int a = 0; a < n_ages; a++) 
      {
        FAA(y,f,a) = F(y,f) * selblocks(selblock_pointer_fleets(y,f)-1,a);
        log_FAA(y,f,a) = log(FAA(y,f,a));
        FAA_tot(y,a) = FAA_tot(y,a) + FAA(y,f,a);
      }
    }
  }

  ZAA = FAA_tot + MAA;
  
  SSB.setZero();
  //year 1
  for(int a = 0; a < n_ages; a++) 
  {
    if(N1_model == 0) NAA(0,a) = exp(log_N1_pars(a));
    else
    {
      if(a==0) NAA(0,0) = exp(log_N1_pars(0));
      else 
      {
        if(a == n_ages-1) NAA(0,a) = NAA(0,a-1)/(1.0 + exp(-MAA(0,a) - exp(log_N1_pars(1)) * FAA_tot(0,a)/FAA_tot(0,n_ages-1)));
        else NAA(0,a) = NAA(0,a-1)* exp(-MAA(0,a) -  exp(log_N1_pars(1)) * FAA_tot(0,a)/FAA_tot(0,n_ages-1));
      }
    }
    SSB(0) += NAA(0,a) * waa(waa_pointer_ssb-1,0,a) * mature(0,a) * exp(-ZAA(0,a)*fracyr_SSB(0));
    pred_NAA(0,a) = NAA(0,a);
  }
  
  //after year 1
  //get predicted numbers at age
  vector<Type> M(n_ages), mat(n_ages), waassb(n_ages), log_SPR0(n_years_model);
  int nh = n_years_model, na = n_years_model;
  vector<Type> log_SR_a(nh), log_SR_b(nh), SR_h(nh), SR_R0(nh);
    
  for(int y = 0; y < n_years_model; y++) 
  {
    for(int a = 0; a < n_ages; a++) 
    {
      M(a) = MAA(y,a);
      waassb(a) = waa(waa_pointer_ssb-1,y,a);
      mat(a) = mature(y,a);
    }
    log_SPR0(y) = log(get_SPR_0(M, mat, waassb, fracyr_SSB(y)));
  }
  REPORT(log_SPR0);
  ADREPORT(log_SPR0);
  if(recruit_model > 2) //BH or Ricker SR
  {
    if(recruit_model == 3) //BH stock recruit
    {
      if(use_steepness == 1)
      {
        SR_h.fill(0.2 + 0.8/(1+exp(-mean_rec_pars(0)))); //SR_a * SPR0/(4.0 + SR_a*SPR0);
        SR_R0.fill(exp(mean_rec_pars(1))); //(SR_a - 1/SPR0) / SR_b;
        log_SR_a = log(4 * SR_h/(exp(log_SPR0)*(1 - SR_h))); //yearly values because of log_SPR0
        if(R0_is_S0 == 1) log_SR_b = log((5*SR_h - 1)/((1-SR_h)*SR_R0)); //yearly values because of log_SPR0
        else log_SR_b = log((5*SR_h - 1)/((1-SR_h)*SR_R0*exp(log_SPR0))); //yearly values because of log_SPR0
      }
      else
      {
        log_SR_a.fill(mean_rec_pars(0));
        log_SR_b.fill(mean_rec_pars(1));
        SR_h = exp(log_SR_a) * exp(log_SPR0)/(4.0 + exp(log_SR_a + log_SPR0));
        if(R0_is_S0 == 1) SR_R0 = (exp(log_SR_a + log_SPR0) - 1) / exp(log_SR_b);
        else SR_R0 = (exp(log_SR_a) - 1/exp(log_SPR0)) / exp(log_SR_b);
      }
    }
    if(recruit_model>3) //Ricker stock recruit
    {
      if(use_steepness == 1)
      {
        SR_h.fill(0.2 + exp(mean_rec_pars(0)));
        SR_R0.fill(exp(mean_rec_pars(1))); 
        log_SR_a = 1.25*log(5*SR_h) - log_SPR0;
        if(R0_is_S0 == 1) log_SR_b = log(1.25*log(5*SR_h)/(SR_R0));
        else log_SR_b = log(1.25*log(5*SR_h)/(SR_R0*exp(log_SPR0)));
      }
      else
      {
        log_SR_a.fill(mean_rec_pars(0));
        log_SR_b.fill(mean_rec_pars(1));
        SR_h = 0.2 * exp(0.8*log(exp(log_SR_a) * exp(log_SPR0)));
        if(R0_is_S0 == 1) SR_R0 = log(exp(log_SR_a + log_SPR0))/(exp(log_SR_b));
        else SR_R0 = log(exp(log_SR_a + log_SPR0))/(exp(log_SR_b + log_SPR0));
      }
    }
    ADREPORT(log_SR_a);
    ADREPORT(log_SR_b);
    vector<Type> logit_SR_h = log(SR_h - 0.2) - log(1 - SR_h);
    vector<Type> log_SR_R0 = log(SR_R0);
    ADREPORT(logit_SR_h);
    ADREPORT(log_SR_R0);
    REPORT(log_SR_a);
    REPORT(log_SR_b);
    REPORT(logit_SR_h);
    REPORT(log_SR_R0);
  }
  matrix<Type> nll_NAA(n_years_model-1,n_ages);
  nll_NAA.setZero();
  vector<Type> nll_recruit(n_years_model-1);
  nll_recruit.setZero();

  for(int y = 1; y < n_years_model; y++) 
  {
    //expected recruitment
    if(recruit_model == 1) pred_NAA(y,0) = NAA(y-1,0); //random walkNAA(y,1)
    else
    {
      if(recruit_model == 2) pred_NAA(y,0) = exp(mean_rec_pars(0)); //random about mean
      else //BH stock recruit
      {
        if(recruit_model == 3) //BH stock recruit
        {
          if(use_steepness == 1) pred_NAA(y,0) = exp(log_SR_a(steepness_year-1)) * SSB(y-1)/(1 + exp(log_SR_b(steepness_year-1))*SSB(y-1));
          else pred_NAA(y,0) = exp(log_SR_a(0)) * SSB(y-1)/(1 + exp(log_SR_b(0))*SSB(y-1));
        }
        else //Ricker stock recruit
        {
          if(use_steepness == 1) pred_NAA(y,0) = exp(log_SR_a(steepness_year-1)) * SSB(y-1) * exp(-exp(log_SR_b(steepness_year-1)) * SSB(y-1)); 
          else pred_NAA(y,0) = exp(log_SR_a(0)) * SSB(y-1) * exp(-exp(log_SR_b(0)) * SSB(y-1));
        }
      }
    }
    //expected numbers at age after recruitment
    for(int a = 1; a < n_ages-1; a++) pred_NAA(y,a) = NAA(y-1,a-1) * exp(-ZAA(y-1,a-1));
    pred_NAA(y,n_ages-1) = NAA(y-1,n_ages-2) * exp(-ZAA(y-1,n_ages-2)) + NAA(y-1,n_ages-1) * exp(-ZAA(y-1,n_ages-1));

    if(use_NAA_re == 1) //random effects NAA, state-space model for all numbers at age
    {
      for(int a = 0; a < n_ages; a++) 
      {
        Type mu = log(pred_NAA(y,a));
        if(bias_correct_pe == 1) mu -= 0.5*exp(2*log_NAA_sigma(NAA_sigma_pointers(a)-1));
        nll_NAA(y-1,a) -= dnorm(log_NAA(y-1,a), mu, exp(log_NAA_sigma(NAA_sigma_pointers(a)-1)), 1);
        SIMULATE 
        {
          if(simulate_state == 1) log_NAA(y-1,a) = rnorm(mu, exp(log_NAA_sigma(NAA_sigma_pointers(a)-1)));
        }
        NAA(y,a) = exp(log_NAA(y-1,a));
      }
    }
    else //recruitments are still estimated parameters, fixed or random effects
    {
      if(random_recruitment == 1) //estimate recruitment as random effects, otherwise as fixed effects. pred_NAA(y,0) should be properly specified above in any case.
      {
        Type mu = log(pred_NAA(y,0));
        if(bias_correct_pe == 1) mu -= 0.5*exp(2*log_R_sigma);
        nll_recruit(y-1) -= dnorm(log_R(y-1), mu, exp(log_R_sigma), 1); 
        SIMULATE if(simulate_state == 1) log_R(y-1) = rnorm(mu, exp(log_R_sigma));
      }
      NAA(y,0) = exp(log_R(y-1));
      //when random effects not used for all numbers at age, survival is deterministic.
      for(int a = 1; a < n_ages; a++) NAA(y,a) = pred_NAA(y,a);
    }
    for(int a = 0; a < n_ages; a++) SSB(y) += NAA(y,a) * waa(waa_pointer_ssb-1,y,a) * mature(y,a) * exp(-ZAA(y,a)*fracyr_SSB(y));
  }
  if(use_NAA_re == 1) 
  {
    REPORT(nll_NAA);
    nll += nll_NAA.sum();
    SIMULATE REPORT(log_NAA);
  }
    
  if(random_recruitment == 1) 
  {
    REPORT(nll_recruit);
    nll += nll_recruit.sum();
    SIMULATE REPORT(log_R);
  }
  
  if(use_b_prior == 1)
  {
    Type mu = log(0.305);
    if(bias_correct_pe == 1) mu -= 0.5*exp(2*log(0.08));
    Type lprior_b = dnorm(log_b, mu, Type(0.08), 1);
    SIMULATE 
    {
      if(simulate_state == 1) log_b = rnorm(mu, Type(0.08));
      REPORT(log_b);
    }
    //see(lprior_b);
    REPORT(lprior_b);
    nll -= lprior_b;
  }
  //see(nll);

  matrix<Type> nll_agg_catch(n_years_catch,n_fleets), nll_catch_acomp(n_years_catch,n_fleets);
  nll_agg_catch.setZero();
  nll_catch_acomp.setZero();
  for(int y = 0; y < n_years_catch; y++)
  {
    int acomp_par_count = 0;
    for(int f = 0; f < n_fleets; f++)
    {
      pred_catch(y,f) = 0.0;
      Type tsum = 0.0;
      for(int a = 0; a < n_ages; a++) 
      {
        pred_CAA(y,f,a) =  NAA(y,a) * FAA(y,f,a) * (1 - exp(-ZAA(y,a)))/ZAA(y,a);
        pred_catch(y,f) += waa(waa_pointer_fleets(f)-1,y,a) * pred_CAA(y,f,a);
        tsum += pred_CAA(y,f,a);
      }
      Type mu = log(pred_catch(y,f));
      Type sig = agg_catch_sigma(y,f)*exp(log_catch_sig_scale(f));
      if(bias_correct_oe == 1) mu -= 0.5*exp(2*log(sig));
      nll_agg_catch(y,f) -= dnorm(log(agg_catch(y,f)), mu, sig,1);
      SIMULATE agg_catch(y,f) = exp(rnorm(mu, sig));
      log_pred_catch(y,f) = log(pred_catch(y,f));
      if(any_fleet_age_comp(f) == 1)
      {
        vector<Type> acomp_pars(n_age_comp_pars_fleets(f));
        for(int j = 0; j < n_age_comp_pars_fleets(f); j++) 
        {
          acomp_pars(j) = catch_paa_pars(acomp_par_count);
          acomp_par_count++;
        }
        if(use_catch_paa(y,f) == 1) 
        {
          for(int a = 0; a < n_ages; a++)
          {
            pred_catch_paa(y,f,a) = pred_CAA(y,f,a)/tsum;
            t_pred_paa(a) = pred_catch_paa(y,f,a);
            t_paa(a) = catch_paa(f,y,a);
          }
          nll_catch_acomp(y,f) -= get_acomp_ll(y, n_ages, catch_Neff(y,f), age_comp_model_fleets(f), t_paa, t_pred_paa, acomp_pars, catch_aref(y,f));
          SIMULATE
          {
            t_paa = sim_acomp(y, n_ages, catch_Neff(y,f), age_comp_model_fleets(f), t_paa, t_pred_paa, acomp_pars, catch_aref(y,f));
            for(int a = 0; a < n_ages; a++) catch_paa(f,y,a) = t_paa(a);
          }
        }
      }
    }
  }
  SIMULATE REPORT(agg_catch);
  SIMULATE REPORT(catch_paa);
  REPORT(nll_agg_catch);
  nll += nll_agg_catch.sum();
  //see(nll);
  REPORT(nll_catch_acomp);
  nll += nll_catch_acomp.sum();
  //see(nll);
  
  matrix<Type> nll_agg_indices(n_years_catch,n_indices), nll_index_acomp(n_years_catch,n_indices);
  nll_agg_indices.setZero();
  nll_index_acomp.setZero();
  pred_indices.setZero();
  for(int y = 0; y < n_years_indices; y++)
  {
    int acomp_par_count = 0;
    for(int i = 0; i < n_indices; i++) 
    {
      Type tsum = 0.0;
      for(int a = 0; a < n_ages; a++) 
      {
        pred_IAA(y,i,a) =  NAA(y,a) * QAA(y,i,a) * exp(-ZAA(y,a) * fracyr_indices(y,i));
        if(units_indices(i) == 1) pred_indices(y,i) += waa(waa_pointer_indices(i)-1,y,a) * pred_IAA(y,i,a);
        else pred_indices(y,i) += pred_IAA(y,i,a);
      }
      for(int a = 0; a < n_ages; a++) 
      {
        if(units_index_paa(i) == 1) pred_IAA(y,i,a) = waa(waa_pointer_indices(i)-1,y,a) * pred_IAA(y,i,a);
        tsum += pred_IAA(y,i,a);
      }
      
      if(use_indices(y,i) == 1)
      {
        Type mu = log(pred_indices(y,i));
        Type sig = agg_index_sigma(y,i)*exp(log_index_sig_scale(i));
        if(bias_correct_oe == 1) mu -= 0.5*exp(2*log(sig));
        nll_agg_indices(y,i) -= dnorm(log(agg_indices(y,i)), mu, sig, 1);
        SIMULATE agg_indices(y,i) = exp(rnorm(mu, sig));
      }
      if(any_index_age_comp(i) == 1)
      {
        vector<Type> acomp_pars(n_age_comp_pars_indices(i));
        for(int j = 0; j < n_age_comp_pars_indices(i); j++) 
        {
          acomp_pars(j) = index_paa_pars(acomp_par_count);
          acomp_par_count++;
        }
        if(use_index_paa(y,i) > 0)
        {
          for(int a = 0; a < n_ages; a++)
          {
            pred_index_paa(y,i,a) = pred_IAA(y,i,a)/tsum;
            t_pred_paa(a) = pred_index_paa(y,i,a);
            t_paa(a) = index_paa(i, y, a);
          }
          nll_index_acomp(y,i) -= get_acomp_ll(y, n_ages, index_Neff(y,i), age_comp_model_indices(i), t_paa, t_pred_paa, acomp_pars, index_aref(y,i));
          SIMULATE
          {
            t_paa = sim_acomp(y, n_ages, index_Neff(y,i), age_comp_model_indices(i), t_paa, t_pred_paa, acomp_pars, index_aref(y,i));
            for(int a = 0; a < n_ages; a++) index_paa(i,y,a) = t_paa(a);
          }
        }
      }
    }
  }
  SIMULATE REPORT(agg_indices);
  SIMULATE REPORT(index_paa);
  REPORT(nll_agg_indices);
  nll += nll_agg_indices.sum();
  //see(nll);
  REPORT(nll_index_acomp);
  nll += nll_index_acomp.sum();
  //see(nll);
  
  
  //////////////////////////////////////////
  //Still need to add in yearly vectors of biological inputs, make sure to calculate SR_a,SR_b vector or otherwise.
  //////////////////////////////////////////
  //calculate BRPs
  //First SPR-based proxies
  //Type percentSPR = 40;
  vector<Type> predR = pred_NAA.col(0);
  matrix<Type> SPR_res = get_SPR_res(MAA, FAA_tot, which_F_age, waa, waa_pointer_ssb, waa_pointer_totcatch, mature, percentSPR, predR, fracyr_SSB, log_SPR0);
  vector<Type> log_FXSPR = SPR_res.col(0);
  vector<Type> log_SSB_FXSPR = SPR_res.col(1);
  vector<Type> log_Y_FXSPR = SPR_res.col(2);
  vector<Type> log_SPR_FXSPR = SPR_res.col(3);
  vector<Type> log_YPR_FXSPR = SPR_res.col(4);
  matrix<Type> log_FXSPR_iter = SPR_res.block(0,5,n_years_model,10);
  REPORT(log_FXSPR_iter);
  REPORT(log_FXSPR);
  REPORT(log_SSB_FXSPR);
  REPORT(log_Y_FXSPR);
  REPORT(log_SPR_FXSPR);
  REPORT(log_YPR_FXSPR);
  ADREPORT(log_FXSPR);
  ADREPORT(log_SSB_FXSPR);
  ADREPORT(log_Y_FXSPR);
  
  //If stock-recruit models
  if(recruit_model > 2) //Beverton-Holt or Ricker selected
  {
    int n = 10;
    vector<Type> log_FMSY(n_years_model), log_FMSY_i(1), waacatch(n_ages), sel(n_ages);
    matrix<Type> log_FMSY_iter(n_years_model,n);
    log_FMSY_iter.col(0).fill(log(FMSY_startval)); //starting value
    vector<Type> log_YPR_MSY(n_years_model), log_SPR_MSY(n_years_model), log_R_MSY(n_years_model);
    Type SR_a, SR_b;
    for(int y = 0; y < n_years_model; y++)
    {
      for(int a = 0; a < n_ages; a++) 
      {
        M(a) = MAA(y,a);
        sel(a) = FAA_tot(y,a)/FAA_tot(y,which_F_age-1); //have to look at FAA_tot to see where max F is.
        waassb(a) = waa(waa_pointer_ssb-1,y,a);
        waacatch(a) = waa(waa_pointer_totcatch-1, y, a);
        mat(a) = mature(y,a);
      }
      if(use_steepness == 1)
      {
        SR_a = exp(log_SR_a(steepness_year-1));
        SR_b = exp(log_SR_b(steepness_year-1));
      }
      else
      {
        SR_a = exp(log_SR_a(0));
        SR_b = exp(log_SR_b(0));
      }
      if(recruit_model == 3) //Beverton-Holt selected 
      {
        sr_yield<Type> sryield(SR_a, SR_b, M, sel, mat, waassb, waacatch,fracyr_SSB(y),0);
        for (int i=0; i<n-1; i++)
        {
          log_FMSY_i(0) = log_FMSY_iter(y,i);
          vector<Type> grad_sr_yield = autodiff::gradient(sryield,log_FMSY_i);
          matrix<Type> hess_sr_yield = autodiff::hessian(sryield,log_FMSY_i);
          log_FMSY_iter(y,i+1) = log_FMSY_iter(y,i) - grad_sr_yield(0)/hess_sr_yield(0,0); 
        }
      }
      else //Ricker selected
      {
        sr_yield<Type> sryield(SR_a, SR_b, M, sel, mat, waassb, waacatch,fracyr_SSB(y),1);
        for (int i=0; i<n-1; i++)
        {
          log_FMSY_i(0) = log_FMSY_iter(y,i);
          vector<Type> grad_sr_yield = autodiff::gradient(sryield,log_FMSY_i);
          matrix<Type> hess_sr_yield = autodiff::hessian(sryield,log_FMSY_i);
          log_FMSY_iter(y,i+1) = log_FMSY_iter(y,i) - grad_sr_yield(0)/hess_sr_yield(0,0); 
        }
      }
      log_FMSY(y) = log_FMSY_iter(y,n-1);
      log_SPR_MSY(y) = log(get_SPR(log_FMSY(y), M, sel, mat, waassb, fracyr_SSB(y)));
      log_YPR_MSY(y) = log(get_YPR(log_FMSY(y), M, sel, waacatch));
      
      if(recruit_model == 3) log_R_MSY(y) = log((SR_a - 1/exp(log_SPR_MSY(y))) / SR_b); //bh
      else log_R_MSY(y) = log(log(SR_a) + log_SPR_MSY(y)) - log(SR_b) - log_SPR_MSY(y); //ricker
    }
    
    vector<Type> log_SSB_MSY = log_R_MSY + log_SPR_MSY;
    vector<Type> log_MSY = log_R_MSY + log_YPR_MSY;
    
    ADREPORT(log_FMSY);
    ADREPORT(log_SSB_MSY);
    ADREPORT(log_R_MSY);
    ADREPORT(log_MSY);
    ADREPORT(log_SPR_MSY);
    ADREPORT(log_YPR_MSY);
    REPORT(log_FMSY);
    REPORT(log_FMSY_iter);
    REPORT(log_SSB_MSY);
    REPORT(log_R_MSY);
    REPORT(log_MSY);
    REPORT(log_SPR_MSY);
    REPORT(log_YPR_MSY);
  }

  matrix<Type> log_FAA_tot = log(FAA_tot.array());
  matrix<Type> log_index_resid = log(agg_indices.block(0,0,n_years_model,n_indices).array()) - log(pred_indices.array());
  matrix<Type> log_catch_resid = log(agg_catch.block(0,0,n_years_model,n_fleets).array()) - log(pred_catch.array());
  vector<Type> log_SSB =  log(SSB);
  vector<Type> Fbar(n_years_model);
  Fbar.setZero();
  int n_Fbar_ages = Fbar_ages.size();
  for(int y = 0; y < n_years_model; y++) for(int a = 0; a < n_Fbar_ages; a++) Fbar(y) += FAA_tot(y,Fbar_ages(a)-1)/n_Fbar_ages;
  
  vector<Type> log_Fbar = log(Fbar);
  matrix<Type> log_NAA_rep = log(NAA.array());

  REPORT(NAA);
  REPORT(pred_NAA);
  REPORT(SSB);
  REPORT(selblocks);
  REPORT(MAA);
  REPORT(q);
  REPORT(QAA);
  REPORT(F);
  REPORT(FAA);
  REPORT(FAA_tot);
  REPORT(Fbar);
  REPORT(pred_catch);
  REPORT(pred_catch_paa);
  REPORT(pred_CAA);
  REPORT(pred_indices);
  REPORT(pred_index_paa);
  REPORT(pred_IAA);
  
  ADREPORT(log_F);
  ADREPORT(log_FAA);
  ADREPORT(log_FAA_tot);
  ADREPORT(log_Fbar);
  ADREPORT(log_NAA_rep);
  ADREPORT(log_SSB);
  ADREPORT(log_pred_catch);
  ADREPORT(pred_IAA);
  ADREPORT(log_index_resid);
  ADREPORT(log_catch_resid);
  
  REPORT(nll);
  
  return nll;
}


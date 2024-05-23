data {
  //dimensions
  int<lower=0> N;          // number of observations
  int<lower=0> p;          // number of predictors
  int<lower=0> j;          // number of census tracts
  //data
  array[N] int<lower=0> y;        // target
  matrix[N, p] X; // matrix of all features
  //Indices
  int<lower=1, upper=j> cusecs_indx[N];  //index of census tract  
  int betas_to_estimate; //number of betas to estimate
  int<lower=1, upper=p> parameter_indx[betas_to_estimate]; //list of the indices of parameters to estimate
  int sig_beta0;  //whether intercept is to be included
}

parameters {
  real beta0; // Intercept
  vector[betas_to_estimate] betas; // Other coefficients

}

model {
  
// Likelihood
for (n in 1:N) {
  
  //1) Calculate the fixed effect 
  real mu_n = 0;
     //Add intercept if needed
  if (sig_beta0==1){
    mu_n += beta0;
  }
     //Go by significant parameter index and add the vector product
  for (par_idx in 1:betas_to_estimate){
    mu_n += X[n,parameter_indx[par_idx]] * betas[par_idx];
  }
  
  //2) Likelihood with randomeffects
  y[n] ~ poisson_log(mu_n);
}


// Wide prior on the significant beats
for (par_idx in 1:betas_to_estimate) {
    betas[par_idx] ~ normal(0, 5);
}

beta0 ~ normal(0, 5);
}


generated quantities {
  // Loglikelihoods and Predictions
  vector[N] log_lik_fixed;
  vector[N] y_pred_fixed_effects;
  for (n in 1:N) {
    real mu_n = 0;
    if (sig_beta0==1){
      mu_n += beta0;
    }
    for (par_idx in 1:betas_to_estimate) {
      mu_n += X[n, parameter_indx[par_idx]] * betas[par_idx];
    }
    
    log_lik_fixed[n] = poisson_log_lpmf(y[n] | mu_n);
    y_pred_fixed_effects[n] = poisson_log_rng(mu_n);
  }
}



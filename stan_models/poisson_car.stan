data {
  //dimensions
  int<lower=0> N;          // number of observations
  int<lower=0> p;          // number of predictors
  int<lower=0> j;          // number of census tracts
  //data
  array[N] int<lower=0> y;        // target
  matrix[N, p] X; // matrix of all features
  //CAR
  int<lower=0> N_edges;
  array[N_edges] int<lower=1, upper=N> node1; // node1[i] adjacent to node2[i]
  array[N_edges] int<lower=1, upper=N> node2; // and node1[i] < node2[i]
  //Indices
  int<lower=1, upper=j> cusecs_indx[N];  //index of census tract  
  int betas_to_estimate; //number of betas to estimate
  int<lower=1, upper=p> parameter_indx[betas_to_estimate]; //list of the indices of parameters to estimate
  int sig_beta0;  //whether intercept is to be included
}

parameters {
  real beta0; // Intercept
  vector[betas_to_estimate] betas; // Other coefficients
  vector[j] phi; // Structured (spatial) effect
  vector[j] theta; // Unstructured (heterogenous) effect
  real<lower=0> tau_theta; // Precision of heterogeneous effects
  real<lower=0> tau_phi; // Precision of spatial effects
  real<lower=0> lambda[p]; // Local shrinkage parameters

}

transformed parameters {
  real<lower=0> sigma_theta = inv(sqrt(tau_theta)); // convert precision to sigma 
    real<lower=0> sigma_phi = inv(sqrt(tau_phi)); // convert precision to sigma
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
  y[n] ~ poisson_log(mu_n + phi[cusecs_indx[n]] * sigma_phi
                + theta[cusecs_indx[n]] * sigma_theta);
}

// Soft sum-to-zero constraint on Phi
target += -0.5 * dot_self(phi[node1] - phi[node2]);
sum(phi) ~ normal(0, 0.001 * j); // equivalent to mean(phi) ~ normal(0,0.001)

// Prior for local shrinkage parameters
lambda ~ cauchy(0, 1);
  
// Wide prior on the significant beats
for (par_idx in 1:betas_to_estimate) {
    betas[par_idx] ~ normal(0, 5);
}

beta0 ~ normal(0, 5);
theta ~ normal(0, 1);
tau_theta ~ gamma(3.2761, 1.81); // Carlin WinBUGS priors
tau_phi ~ gamma(1, 1); // Carlin WinBUGS priors
}


generated quantities {
  // Loglikelihoods and Predictions
  vector[N] log_lik;
  vector[N] log_lik_fixed;
  vector[N] y_pred;
  vector[N] y_pred_fixed_and_icar;
  vector[N] y_pred_fixed_effects;
  for (n in 1:N) {
    real mu_n = 0;
    if (sig_beta0==1){
      mu_n += beta0;
    }
    for (par_idx in 1:betas_to_estimate) {
      mu_n += X[n, parameter_indx[par_idx]] * betas[par_idx];
    }
    real spatial_effect = phi[cusecs_indx[n]] * sigma_phi;
    real heterogeneous_effect = theta[cusecs_indx[n]] * sigma_theta;
    
    log_lik[n] = poisson_log_lpmf(y[n] | mu_n + spatial_effect + heterogeneous_effect);
    log_lik_fixed[n] = poisson_log_lpmf(y[n] | mu_n);
    y_pred[n] = poisson_log_rng(mu_n + spatial_effect + heterogeneous_effect);
    y_pred_fixed_and_icar[n] = poisson_log_rng(mu_n + spatial_effect);
    y_pred_fixed_effects[n] = poisson_log_rng(mu_n);
  }
}



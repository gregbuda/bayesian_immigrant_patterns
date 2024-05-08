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
}

parameters {
  real beta0; // Intercept
  vector[j] phi; // Structured (spatial) effect
  vector[j] theta; // Unstructured (heterogenous) effect
  real<lower=0> tau_theta; // Precision of heterogeneous effects
  real<lower=0> tau_phi; // Precision of spatial effects
}

transformed parameters {
  real<lower=0> sigma_theta = inv(sqrt(tau_theta)); // convert precision to sigma 
    real<lower=0> sigma_phi = inv(sqrt(tau_phi)); // convert precision to sigma 
}


model {
  // Likelihood
  for (n in 1:N) {
    y[n] ~ poisson_log(beta0  + phi[cusecs_indx[n]] * sigma_phi
                  + theta[cusecs_indx[n]] * sigma_theta);
  }
  
  // Soft sum-to-zero constraint on Phi
    target += -0.5 * dot_self(phi[node1] - phi[node2]);
    sum(phi) ~ normal(0, 0.001 * j); // equivalent to mean(phi) ~ normal(0,0.001)
  
  // Priors
  beta0 ~ normal(0, 5);
  theta ~ normal(0, 1);
  tau_theta ~ gamma(3.2761, 1.81); // Carlin WinBUGS priors
  tau_phi ~ gamma(1, 1); // Carlin WinBUGS priors
}

generated quantities {
  // 1a) Likelihoods
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = poisson_log_lpmf(y[n] | beta0 + phi[cusecs_indx[n]] * sigma_phi
                                + theta[cusecs_indx[n]] * sigma_theta);
  }
  
    // 1b) Likelihoods only for fixed effects
  vector[N] log_lik_fixed;
  for (n in 1:N) {
    log_lik_fixed[n] = poisson_log_lpmf(y[n] | beta0);
  }

  // 2) Predictions for y
  int y_pred[N];
  for (n in 1:N) {
    y_pred[n] = poisson_log_rng(beta0 + phi[cusecs_indx[n]] * sigma_phi
                                + theta[cusecs_indx[n]] * sigma_theta);
  }
  
    // 3) Predictions for y without unstructured
  int y_pred_icar[N];
  for (n in 1:N) {
    y_pred_icar[n] = poisson_log_rng(beta0 + phi[cusecs_indx[n]] * sigma_phi);
  }
}
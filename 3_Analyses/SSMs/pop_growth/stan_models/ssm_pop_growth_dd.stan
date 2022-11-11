//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

functions {
}

data {
  int Time;
  vector<lower = 0>[Time] logy;
  real sigma_obs;
  vector<lower = 0>[Time] prec;
  vector[Time] tmin;
  vector[Time] tmax;
  int<lower = 0, upper = 1> dry[Time];
}

transformed data {
  real<lower = 0> logy1 = logy[1];
}

parameters {
  real logN_est1; // Initial log population size
  real b; // Mean growth rate
  real dd; // Density-dependent coefficient
  real<lower = 0, upper = 1> sigma_proc; // SD of state procesdatas
  // real<lower = 0, upper = 1> sigma_obs; // SD of observation process
  vector[Time - 1] r;
}

transformed parameters {
  vector[Time] logN_est; // Log population size
  
  // State process
  logN_est[1] = logN_est1;
  for (t in 1 : (Time - 1)) {
    logN_est[t + 1] = logN_est[t] - dd*exp(logN_est[t]) + r[t];
  }
}

model {
  // Priors and constraints
  logN_est1 ~ normal(logy1, sqrt(10));
  b ~ normal(1, sqrt(1000));
  dd ~ gamma(0.25, 0.25);
  //  sigma_proc ~ uniform(0, 1);
  //  sigma_obs ~ uniform(0, 1);
  
  // Likelihood
  r ~ normal(b, sigma_proc);
  
  // Observation process
  logy ~ normal(logN_est, sigma_obs);
}

generated quantities {
  real sigma2_proc;
  // real sigma2_obs;
  vector<lower = 0>[Time] N_est; // Population size
  
  sigma2_proc = square(sigma_proc);
  // sigma2_obs = square(sigma_obs);
  
  for (t in 1:Time) {
    N_est[t] = exp(logN_est[t]);
  }
  
}

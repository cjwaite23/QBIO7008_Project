// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// Jolly-Seber model as a restricted occupancy model
// code modified from Ito, 2022

functions {
  // These functions are derived from Section 12.3 of
  // Stan Modeling Language User's Guide and Reference Manual
  
  /**
   * Return a integer value of first capture occasion
   *
   * @param y_i Integer array of capture history
   *
   * @return Integer value of first capture occasion
   */
  int first_cap(int[] y_i) {
    if (y_i[sort_indices_desc(y_i)[1]]) {
      return sort_indices_desc(y_i)[1];
      } else {
        return 0;
      }
    }
  
  /**
   * Return a integer value of last capture occasion
   *
   * @param y_i Integer array of capture history
   *
   * @return Integer value of last capture occasion
   */
  int last_cap(int[] y_i, int n_occasions) {
    if (y_i[sort_indices_asc(y_i)[n_occasions]]) {
      return sort_indices_asc(y_i)[n_occasions];
      } else {
        return 0;
      }
    }
  
  /**
   * Return a matrix of uncaptured probabilities
   *
   * @param p           Matrix of detection probabilities for each individual
   *                    and capture occasion
   * @param phi         Matrix of survival probabilities for each individual
   *                    and capture occasion
   *
   * @return Vector of uncaptured probabilities
   */
  vector prob_uncaptured(vector p, vector phi) {
    int n_occasions = rows(p);
    vector[n_occasions] chi;
    
    chi[n_occasions] = 1.0;
    for (t in 1 : (n_occasions - 1)) {
      int t_curr = n_occasions - t;
      int t_next = t_curr + 1;
        
      chi[t_curr] = (1 - phi[t_curr]) + phi[t_curr] * (1 - p[t_next])
                      * chi[t_next];
    }
    return chi;
  }
  
  /**
   * Calculate log likelihood of a Jolly-Seber model
   *
   * @param y     Integer array of capture history
   * @param first Integer array of first capture occasions
   * @param last  Integer array of last capture occasions
   * @param p     Vector of detection probabilities
   * @param phi   Vector of survival probabilities
   * @param gamma Vector of removal entry probabilities
   * @param chi   Vector of uncapture probabilities
   */
  void jolly_seber_lp(int[,] y, int[] first, int[] last,
                      vector p, vector phi, vector gamma, vector chi) {
    int n_ind = dims(y)[1];
    int n_occasions = dims(y)[2];
    vector[n_occasions] qgamma = 1.0 - gamma;
    vector[n_occasions] qp = 1.0 - p;
    
    for (i in 1 : n_ind) {
      if (first[i] == 1) {
      1 ~ bernoulli(gamma[1] * p[1]);
      // Until last capture
      for (t in (first[i] + 1) : last[i]) {
        1 ~ bernoulli(phi[t - 1]); // Survived
        y[i, t] ~ bernoulli(p[t]); // Capture/Non-capture
      }
      // Subsequent occasions
      1 ~ bernoulli(chi[last[i]]);
    } else if (first[i] > 1) {
      vector[first[i]] lp;
          
      // Entered at 1st occasion
      lp[1] = log(gamma[1]) +  
              log(prod(qp[1 : (first[i] - 1)])) +
              log(prod(phi[1 : (first[i] - 1)])) +
              log(p[first[i]]);
      // Entered at t-th occasion (1 < t < first[i])
      for (t in 2 : (first[i] - 1)) {
        lp[t] = log(prod(qgamma[1 : (t - 1)])) +
                log(gamma[t]) +
                log(prod(qp[t : (first[i] - 1)])) +
                log(prod(phi[t : (first[i] - 1)])) +
                log(p[first[i]]);
      }
      lp[first[i]] = log(prod(qgamma[1 : (first[i] - 1)])) +
                     log(gamma[first[i]]) +
                     log(p[first[i]]);
      target += log_sum_exp(lp);
          
      // Until last capture
      for (t in (first[i] + 1) : last[i]) {
        1 ~ bernoulli(phi[t - 1]); // Survived
        y[i, t] ~ bernoulli(p[t]); // Capture/Non-capture
      }
      // Subsequent occasions
      1 ~ bernoulli(chi[last[i]]);
      } else {
        // Never captured
        vector[n_occasions + 1] lp;
        
        // Entered at 1st occasion, but never captured
        lp[1] = log(gamma[1] * qp[1] * chi[1]);
        // Entered at t-th occasion, but never captured
        for (t in 2 : n_occasions) {
          lp[t] = log(prod(qgamma[1 : (t - 1)]) * gamma[t] * qp[t] * chi[t]);
        }
        // Never entered
        lp[n_occasions + 1] = log(prod(qgamma));
        target += log_sum_exp(lp);
      }
    }
  }
  
  /**
   * Returns delta, where
   * delta[n] = gamma[n] * PROD_{m < n} (1 - gamma[m])
   * (Thanks to Dr. Carpenter)
   *
   * @param gamma Vector of probability sequence
   *
   * @return Vector of complementary probability sequence
   */
  vector seq_cprob(vector gamma) {
    int N = rows(gamma);
    vector[N] log_cprob;
    real log_residual_prob = 0;
    
    for (n in 1 : N) {
      log_cprob[n] = log(gamma[n]) + log_residual_prob;
      log_residual_prob = log_residual_prob + log(1 - gamma[n]);
    }
    return exp(log_cprob);
  }
}
data {
  int<lower=0> M; // Augmented sample size
  int<lower=0> n_occasions; // Number of capture occasions
  int<lower=0, upper=1> y[M, n_occasions]; // Augmented capture-history
}
transformed data {
  int n_occ_minus_1 = n_occasions - 1;
  int<lower=0, upper=n_occasions> first[M];
  int<lower=0, upper=n_occasions> last[M];
  
  for (i in 1 : M) {
    first[i] = first_cap(y[i]);
  }
  for (i in 1 : M) {
    last[i] = last_cap(y[i], n_occasions);
  }
}
parameters {
  vector<lower=0, upper=1>[n_occ_minus_1] phi; // Survival
  vector<lower=0, upper=1>[n_occasions] p; // Capture probability
  vector<lower=0, upper=1>[n_occasions] gamma; // Removal entry probability
}
transformed parameters {
  vector<lower=0, upper=1>[n_occasions] chi;

  chi = prob_uncaptured(p, phi);
}
model {
  // Priors
  // Uniform priors are implicitly defined.
  
  // Likelihood
  jolly_seber_lp(y, first, last, p, phi, gamma, chi);
}
generated quantities {
  real psi; // Inclusion probability
  vector[n_occasions] b; // Entry probability
  int Nsuper; // Superpopulation size
  int N[n_occasions]; // Actual population size
  int B[n_occasions]; // Number of entries
  int z[M, n_occasions]; // Latent state
  
  // Generate z[]
  for (i in 1 : M) {
    int q = 1;
    real mu2;
    
    z[i, 1] = bernoulli_rng(gamma[1]);
    for (t in 2 : n_occasions) {
      q = q * (1 - z[i, t - 1]);
      mu2 = phi[t - 1] * z[i, t - 1] + gamma[t] * q;
      z[i, t] = bernoulli_rng(mu2);
    }
  }
  
  // Calculate derived population parameters
  {
    vector[n_occasions] cprob = seq_cprob(gamma);
    int recruit[M, n_occasions] = rep_array(0, M, n_occasions);
    int Nind[M];
    int Nalive[M];
    
    psi = sum(cprob);
    b = cprob / psi;
    
    for (i in 1 : M) {
      int f = first_cap(z[i,  : ]);
      
      if (f > 0) {
        recruit[i, f] = 1;
      }
    }
    for (t in 1 : n_occasions) {
      N[t] = sum(z[ : , t]);
      B[t] = sum(recruit[ : , t]);
    }
    for (i in 1 : M) {
      Nind[i] = sum(z[i]);
      Nalive[i] = 1 - !Nind[i];
    }
    Nsuper = sum(Nalive);
  }
}

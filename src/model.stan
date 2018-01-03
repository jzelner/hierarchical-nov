data {
  int C; //Number of camps
  int T; //Number of days observed
  int N; //Total number of observations
  vector[C] P; //Population of each camp
  int<lower=1, upper=C> J[N]; //Camp index for each observation
  int<lower=1, upper=T> t[N]; //Day of each observation
  vector<lower=0>[N] Y; //Number of cases observed on each day
}

transformed data {
  vector<lower=0>[T] AY = rep_vector(0, T); //Total number of cases on each day
  matrix<lower=0>[C,T] Ymat = rep_matrix(0, C, T); //Total number of cases in each camp
  vector<lower=0>[C] S; //Number of survivors in each camp at end of outbreak
  real day_delta = 1e-9;
  real day_x[T];
  for (i in 1:N) {
    AY[t[i]] = AY[t[i]] + Y[i];
    Ymat[J[i],t[i]] = Ymat[J[i],t[i]] + Y[i];
  }
  {
    real d = 1;
    for (i in 1:(T)) {
    day_x[i] = d;
    d = d + 1;
    }
  }

  for (i in 1:C) {
    S[i] = P[i] - sum(Ymat[i]);
  }


}

parameters {
  real log_beta_mu; //Avg log beta
  matrix<lower=0>[C,T] beta; //Realized total beta by day
  real<lower=0> beta_shape; //Shape of distribution of betas
  real<lower=0, upper=1> alpha; //Per-capita exposure to individuals outside camp
  real<lower=0> gamma; //Shape of infectious period

}

transformed parameters {
  matrix[C,T] lambda = rep_matrix(0, C, T); //Daily per-capita force of infection for each camp
  matrix[C,T] c_lambda = rep_matrix(0, C, T);

  vector<lower=0, upper=1>[T] inf_day; //Distribution of infectiousness by day
  real beta_rate = beta_shape/exp(log_beta_mu);

  //Pre-calculate proportion of infectiousness on each day since onset
  inf_day[1] = exponential_cdf(1, gamma);
  for (i in 2:T) {
    inf_day[i] = exponential_cdf(i, gamma) - exponential_cdf(i-1, gamma);
  }

  //Sum across camps to get force of infection for each day
  for (c in 1:C) { //Repeat for each camp
    for (tb in 1:T) {
      for (te in tb:T) {
        int dayindex = te-tb+1;
        real b_t = beta[c,tb];
        real a_t = alpha*b_t;
        real incamp = (b_t/P[c])*Ymat[c,tb]*inf_day[dayindex];
        real outcamp = (a_t/sum(P)) * ((AY[tb]-Ymat[c,tb])*inf_day[dayindex]);
        lambda[c,te] = lambda[c, te] + ((incamp + outcamp));
      }
    }
  }

  //Now calculate cumulative force of infection for individuals
  //infected on each day
  for (c in 1:C) {
    c_lambda[c] = cumulative_sum(lambda[c]);
  }
 
}

model {
  //Model pars
  to_vector(beta) ~ gamma(beta_shape, beta_rate);
  log_beta_mu ~ normal(0, 1);
  gamma ~ normal(0,1);

  //Iterate over camps
  for (c in 1:C) {
    //Log-likelihood for survival
    target += -S[c]*c_lambda[c,T];
    for (i in 2:T) {
      //Log-likelihood for infections
      real c_ex = 0;
      real d_ex = lambda[c,i-1];
      real ll;
      if (i > 2) {
        c_ex = c_lambda[c,i-2];
      }

      if (i > 3) {
        real c_ex_1 = c_lambda[c, i-3];
        real d_ex_1 = lambda[c, i-2];
        ll = log_sum_exp(log(0.5) + log(d_ex)-c_ex, log(0.5) + log(d_ex_1)-c_ex_1);
      } else {
        ll = log(d_ex)-c_ex;
      }
      target += Ymat[c,i]*ll;
    }
  }

}


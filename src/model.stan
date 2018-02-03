data {
  int C; //Number of camps
  int T; //Number of days observed
  int N; //Total number of observations
  vector[C] P; //Population of each camp
  int<lower=1, upper=C> J[N]; //Camp index for each observation
  int<lower=1, upper=T> t[N]; //Day of each observation
  vector<lower=0>[N] Y; //Number of cases observed on each day
  real<lower=0, upper=1> epsilon; //Rate of progression from infection to symptoms
}

transformed data {
  matrix<lower=0>[C,T] Ymat = rep_matrix(0, C, T); //Total number oflow cases in each camp
  vector<lower=0>[C] S; //Number of survivors in each camp at end of outbreak
  vector[T] e_log_ll;
  real total_pop = sum(P);
  for (i in 1:N) {
    Ymat[J[i],t[i]] = Y[i];
  }


  //Calculate number of survivors at end
  for (i in 1:C) {
    S[i] = P[i] - sum(Ymat[i]);
  }

  //Pre-load latency log-lls
  for (i in 1:T) {
    e_log_ll[i] = log(pow((1-epsilon),i-1)*epsilon);
  }
}

parameters {
  real log_beta_mu; //Avg log beta
  vector<lower=0>[N] beta; //Realized avg beta by day
  real<lower=0> beta_shape; //Shape of distribution of betas
  real<lower=0> zeta; //Per-capita exposure to individuals outside camp
  real<lower=0, upper = 1> gamma; //Shape of infectious period

}

transformed parameters {
  matrix[C,T] lambda = rep_matrix(0, C, T); //Daily per-capita force of infection for each camp
  matrix[C,T] c_lambda = rep_matrix(0, C, T);
  matrix[C,T] beta_mat = rep_matrix(0, C, T);
  vector<lower=0, upper=1>[T] inf_day; //Distribution of infectiousness by day

  for (i in 1:N) {
      beta_mat[J[i],t[i]] = beta[i];
  }
  //Pre-calculate proportion of infectiousness on each day since onset
  inf_day[1] = gamma; //exponential_cdf(1, gamma);
  for (i in 2:T) {
    inf_day[i] = pow((1-gamma),i-1)*gamma;//exponential_cdf(i, gamma) - exponential_cdf(i-1, gamma);
  }

  //Sum across camps to get force of infection for each day
  for (c in 1:C) { //Repeat for each camp
    for (c2 in 1:C) {
      real a = c == c2 ? zeta : 1;
      for (tb in 1:T) {
        real b_t = a*beta_mat[c,tb];
        for (te in tb:T) {
          int dayindex = te-tb+1;
          lambda[c2,te] = lambda[c2,te] + b_t*Ymat[c,tb]*inf_day[dayindex];
        }
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
  beta ~ gamma(beta_shape, beta_shape/exp(log_beta_mu));
  log_beta_mu ~ normal(0, 2);
  
  beta_shape ~ normal(1,1);
  #zeta ~ normal(1, 5);
  //Iterate over camps
  for (c in 1:C) {
    //Log-likelihood for survival
    target += -S[c]*c_lambda[c,T];
  }

  for (i in 1:N) {
    if (t[i] > 1) {
      int max_inf_t = t[i]-1; //> 11 ? 11 : t[i]-1;
      vector[max_inf_t] log_ll;
      vector[max_inf_t] elogs;
      real total_e;
      for (tt in 1:max_inf_t) {//tt = possible infection time
      //Need to change this to matrix slice
        real c_ex = tt > 1 ? c_lambda[J[i], tt-1] : 0;
        real d_ex = 1.0 - exp(-lambda[J[i],tt]);
        elogs[tt] = e_log_ll[t[i]-tt];
        log_ll[tt] = log(d_ex) - c_ex;
      }

      total_e = log_sum_exp(elogs);
      for (tt in 1:max_inf_t) {
        log_ll[tt] = log_ll[tt] + elogs[tt] - total_e;
      }

      target += log_sum_exp(log_ll);
    }
  }

  

}

generated quantities {

  vector[T] daily_avg_r;
  matrix[C,T] camp_r;
  matrix[C,T] beta_Y;
  vector[C] c_weights = P / total_pop;
  {

    for (c in 1:C) {
      for (tt in 1:T) {
        real active_y = Ymat[c,tt] > 0 ? beta_mat[c,tt] : 0;
        camp_r[c,tt] = P[c]*(zeta*active_y) + (total_pop - P[c])*active_y;
      }
    }
    beta_Y = camp_r .* Ymat;
    for (tt in 1:T) {
      daily_avg_r[tt] = sum(col(beta_Y,tt))/sum(col(Ymat,tt));
    }



  }
}

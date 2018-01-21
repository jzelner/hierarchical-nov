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
  vector<lower=0>[T] AY = rep_vector(0, T); //Total number of cases on each day
  matrix<lower=0>[C,T] Ymat = rep_matrix(0, C, T); //Total number of cases in each camp
  vector<lower=0>[C] S; //Number of survivors in each camp at end of outbreak
  matrix[C,C] p_pop;
  vector[T] e_log_ll;
  real total_pop = sum(P);
  for (i in 1:N) {
    AY[t[i]] = AY[t[i]] + Y[i];
    Ymat[J[i],t[i]] = Ymat[J[i],t[i]] + Y[i];
  }


  for (i in 1:C) {
    S[i] = P[i] - sum(Ymat[i]);
  }

  for (i in 1:C) {
    for (j in 1:C) {
      if (i == j) {
        p_pop[i,j] = 1;
      } else {
        p_pop[i,j] = P[j]/(sum(P)-P[i]);
      }
    }
  }

  for (i in 1:T) {
    e_log_ll[i] = log(pow((1-epsilon),i-1)*epsilon);
  }
}

parameters {
  vector[2] log_beta; //Avg log beta
  matrix<lower=0>[C,T] beta; //Realized total beta by day
  real<lower=0> beta_shape; //Shape of distribution of betas
  real<lower=0> zeta; //Per-capita exposure to individuals outside camp
  real<lower=0, upper = 1> gamma; //Shape of infectious period

}

transformed parameters {
  matrix[C,T] lambda = rep_matrix(0, C, T); //Daily per-capita force of infection for each camp
  matrix[C,T] c_lambda = rep_matrix(0, C, T);

  vector<lower=0, upper=1>[T] inf_day; //Distribution of infectiousness by day
  //Pre-calculate proportion of infectiousness on each day since onset
  inf_day[1] = gamma; //exponential_cdf(1, gamma);
  for (i in 2:T) {
    inf_day[i] = pow((1-gamma),i-1)*gamma;//exponential_cdf(i, gamma) - exponential_cdf(i-1, gamma);
  }

  //Sum across camps to get force of infection for each day
  for (c in 1:C) { //Repeat for each camp
    for (c2 in 1:C) {
      real a = c == c2 ? zeta : 1;
      for (tb in 1:11) {
        for (te in tb:11) {
          int dayindex = te-tb+1;
          real b_t = a*beta[c,tb]/total_pop;
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

  vector[T] log_ll;
  //Model pars
  for (i in 1:T) {
    real beta_rate = exp(log_beta[1] + log_beta[2]*i)/beta_shape;
    col(beta,i) ~ gamma(beta_shape, beta_rate);
  }
  log_beta ~ normal(0, 1);
  beta_shape ~ normal(1,1);
  gamma ~ normal(0,1);
  zeta ~ normal(0, 1);
  //Iterate over camps
  for (c in 1:C) {
    //Log-likelihood for survival
    target += -S[c]*c_lambda[c,T];
      for (tt in 2:T) {
        if (Ymat[c,tt] > 0) {
        //Log-likelihood for infections
          for (i in 1:(tt-1)) { //i indexes infection times, tt onset times
          real c_ex = i == 1 ? 0 : c_lambda[c,i-1];
          real d_ex = lambda[c,i];
          real e_log = e_log_ll[tt-i];
          log_ll[i] = log(d_ex) - c_ex + e_log;
        }
          target += Ymat[c,tt]*log_sum_exp(head(log_ll, tt-1));
      }
    }
  }

}

generated quantities {

  vector[T] daily_avg_r;
  matrix[C,T] camp_r;
  matrix[C,T] beta_Y;
  {

    for (c in 1:C) {
      for (tt in 1:T) {
        real active_y = Ymat[c,tt] > 0 ? beta[c,tt] : 0;
        real c_p = (P[c]/total_pop);
        camp_r[c,tt] = (c_p)*(zeta*active_y) + (1-c_p)*active_y;
      }
    }
    beta_Y = camp_r .* Ymat;
    for (tt in 1:T) {
      daily_avg_r[tt] = sum(col(beta_Y,tt))/sum(col(Ymat,tt));
    }



  }
}

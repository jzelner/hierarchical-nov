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
  matrix[C,C] bc_pop;
  for (i in 1:N) {
    Ymat[J[i],t[i]] = Y[i];
  }
  //bc_pop is the contact rate with individuals 
  //outside of the case's camp
  for (i in 1:C) {
    for (j in i:C) {
     bc_pop[i,j] = i == j ? 0 : 1/(sum(P)-P[j]);
     bc_pop[j,i] = j == i ? 0 : 1/(sum(P)-P[i]);
    }
  }

  //Calculate number of survivors at end
  for (i in 1:C) {
    S[i] = P[i] - sum(Ymat[i]);
  }

  //Pre-load latency log-lls
  e_log_ll[1] = lognormal_lcdf(1 | 1, 1.5);
  for (i in 2:T) {
#    e_log_ll[i] = log(pow((1-epsilon),i-1)*epsilon);
    e_log_ll[i] = log(lognormal_cdf(i, 1, 1.5)-lognormal_cdf(i-1, 1, 1.5));
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
  matrix[C,T] lambda_w = rep_matrix(0, C, T);
  matrix[C,T] lambda_b = rep_matrix(0, C, T);
  matrix[C,T] c_lambda = rep_matrix(0, C, T);
  matrix[C,T] beta_mat = rep_matrix(0, C, T);
  vector<lower=0, upper=1>[T] inf_day; //Distribution of infectiousness by day
  matrix[C,T] zeta_t = zeta * (bc_pop * Ymat);

   for (i in 1:N) {
       beta_mat[J[i],t[i]] = beta[i]/P[J[i]];
   }
  //Pre-calculate geometrically-distributed proportion of infectiousness on each day since onset
  inf_day[1] = gamma; 
  for (i in 2:T) {
    inf_day[i] = pow((1-gamma),i-1)*gamma;
  }

  //Sum across camps to get force of infection for each day
  for (c in 1:C) { //Repeat for each camp
      for (tb in 1:T) {
        real b_t = beta_mat[c,tb];
        real db = (b_t*Ymat[c,tb]);
        //Get contributions from outside
        real dz = zeta_t[c,tb];
        
        for (te in tb:T) {
          int dayindex = te-tb+1;
          lambda_w[c,te] = lambda_w[c, te] + db*inf_day[dayindex];
          lambda_b[c,te] = lambda_b[c, te] + dz*inf_day[dayindex];
          lambda[c,te] = lambda[c,te] + lambda_w[c, te] + lambda_b[c, te];
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

  //Prior for total infectiousness for each day/camp combination
  for (i in 1:N) {
    beta[i] ~ gamma(beta_shape*Y[i], (beta_shape*Y[i])/exp(log_beta_mu));
  }
  
  
  log_beta_mu ~ normal(-1, 2);
  beta_shape ~ normal(4,1);
  
  //Iterate over camps
  for (c in 1:C) {
    //Log-likelihood for survival
    target += -S[c]*c_lambda[c,T];
  }

  for (i in 1:N) { //Looping over all days where there are > 0 cases
    if (t[i] > 1) {
      int max_inf_t = t[i]-1; 
      vector[max_inf_t] log_ll;
      vector[max_inf_t] elogs;
      real total_e;
      
      //Pre-calculating the LL contributions from 
      //exposure and latency for each possible 
      //day of infection
      for (tt in 1:max_inf_t) {//tt = possible infection time
        //0  cumulative exposure (c_ex) for cases infected on day 0
        real c_ex = tt > 1 ? c_lambda[J[i], tt-1] : 0;  
        real d_ex = 1.0 - exp(-lambda[J[i],tt]);
        elogs[tt] = e_log_ll[t[i]-tt];
        log_ll[tt] = log(d_ex) - c_ex;
      }

      //Adding in log-likelihood associated with latent period
      total_e = log_sum_exp(elogs);
      for (tt in 1:max_inf_t) {
        log_ll[tt] = log_ll[tt] + elogs[tt]; //- total_e;
      }

      target += Y[i]*log_sum_exp(log_ll);
    }
  }

  

}

generated quantities {

  vector[T] daily_avg_r;
  matrix[C,T] camp_r;
  matrix[C,T] beta_Y;
  matrix[C,T] lambda_within;
  vector[C] c_weights = P / total_pop;
  {

    for (c in 1:C) {
      real real_t = 1;
      for (tt in 1:T) {
        real active_y = Ymat[c,tt] > 0 ? 1 : 0;
        camp_r[c,tt] = active_y*(P[c]*beta_mat[c, tt] + (total_pop - P[c])*(zeta/(sum(P)-P[c])));
        lambda_within[c,tt] = lambda_w[c,tt]/(lambda_b[c,tt]+lambda_w[c,tt]);
      }
    }
    beta_Y = camp_r .* Ymat;
    for (tt in 1:T) {
      daily_avg_r[tt] = sum(col(beta_Y,tt))/sum(col(Ymat,tt));
    }



  }
}

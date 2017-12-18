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
  matrix<lower=0>[C,T] PAR = rep_matrix(0, C, T); //Population at risk in each camp
  for (i in 1:N) {
    AY = AY[t[i]] + Y[i];
    Ymat = Ymat[J[i],t[i]] + Y[i];
  }

  //Get population at risk *at time of infection* for incident cases
  //on day t
  for (c in 1:C) {
    PAR[c,1] = P[c];
    for (t in 2:T) {
      PAR[c,t] = PAR[c,t-1]-AY[c,t-1];
    }
  }

}

parameters {
  real log_mu_beta; //Avg log beta
  vector[C] log_beta; //Realized log betas
  real<lower=0> log_beta_sigma; //Variance of log betas
  real<lower=0> alpha; //Per-capita exposure to individuals outside camp
  real<lower=0> gamma; //Shape of infectious period
}

transformed parameters {
  matrix[C,T] lambda = rep_matrix(0, C, T); //Daily per-capita force of infection for each camp
  matrix[C,T] c_lambda = rep_matrix(0, C, T);

  vector<lower=0, upper=>[T] inf_day; //Distribution of infectiousness by day
  vector<lower=0>[C] beta = exp(log_beta);

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
        real incamp = beta[c]*Ymat[c,tb]*day_inf[dayindex];
        real outcamp = alpha * ((AY[tb]-Ymat[c,tb])*day_inf[dayindex]);
        lambda[c,te] = lambda[c, te] + ((incamp + outcamp) / P[c]);
      }
    }
  }

  //Now calculate cumulative force of infection for individuals
  //infected on each day
}

model {

  //Sample log-betas from normal distribution
  log_beta ~ normal(log_mu_beta, sigma_beta);

  for (c in 1:C) {
    for (t in 2:T) {
      target += 
    }
  }

)


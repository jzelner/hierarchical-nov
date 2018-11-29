data {
  
  int S; ## Number of survivors at end of cruise
  int T; ## Length of outbreak
  vector Y[T]; ## Number of cases on each day
  vector[T] SSE; ## Number of SSEs on each day
  matrix[T,T] e_weights;
  matrix[T,T] inst_exposure;
  matrix[T,T] c_exposure;
  
}

parameters {
  
  real mu_r;
  real k;
  real gamma;
  
}
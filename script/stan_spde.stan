
//========================
// DATA
//========================
data {
  int n;    // n obs
  int p;    // n par
  int row_spar;    // rows sparse matrix
  int col_spar;    // cols sparse matrix
  
  vector[n] time;         //The response
  int notcens[n];         //indicator vector stating which persons that were censored
  matrix[n, p] X;         //Design matrix for fixed effects
  matrix[row_spar, col_spar] M0;     // SPDE matrices from INLA
  matrix[row_spar, col_spar] M1;
  matrix[row_spar, col_spar] M2;
  matrix[n, col_spar] A;     //Matrix for interpolating points witin triangles
}

//========================
// PARAMETERS
//========================
parameters {
  vector[p] beta;  
  real log_tau;
  real log_kappa;
  real log_omega;  
  vector[row_spar] u;   // spatial random effect  
 }


//========================
// T.PARAMETERS
//========================
transformed parameters {
  real tau = exp(log_tau);
  real kappa = exp(log_kappa);
  real omega = exp(log_omega);  // Parameter of Weibull distribution
//------------------------------------------
  vector[n] eta;
  matrix[row_spar, col_spar] Q;
  vector[row_spar] zeroes = rep_vector(0, col_spar);
  vector[n] delta;
  delta = (A*u)/tau;
  eta = X*beta + delta;
  Q = kappa^4*M0 + 2.0*kappa^2*M1 + M2;
   
}


//========================
// MODEL
//========================
model {
  for(i in 1:n){    
    real lambda = exp(eta[i]);
    real t_omega = pow(time[i],omega);
    real S = exp(-lambda*t_omega);            // Survival function
    real f = lambda*omega*t_omega/time[i]*S;  // Weibull density
    if(notcens[i]){
      target += log(f);
    }else{
      target += log(S); // The pasient survived until cencoring
    }
  }
  //---------------------------------------------
  u ~ multi_normal_prec(zeroes, Q); // Spatial random effect with precision matrix Q
  log_tau ~ normal(-1, 0.5);
  log_kappa ~ normal(1, 0.5);
}

//========================
// G.QUANTITIES
//========================
generated quantities{
  real range = sqrt(8)/kappa;
}


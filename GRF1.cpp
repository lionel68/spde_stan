
// include libraries
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
#include <string>
using namespace density;
using Eigen::SparseMatrix;
// Can also choose which likelihood to use.
// Lognormal density
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}
// Inverse gamma
template<class Type>
Type dinvgauss(Type x, Type mean, Type shape, int give_log=0){
  Type logres = 0.5*log(shape) - 0.5*log(2*M_PI*pow(x,3)) - (shape * pow(x-mean,2) / (2*pow(mean,2)*x));
  if(give_log) return logres; else return exp(logres);
}
// dcauchy for hyperparameters
template<class Type>
Type dcauchy(Type x, Type mean, Type shape, int give_log=0){
  Type logres = 0.0;
  logres-= log(M_PI);
  logres-= log(shape);
  // Note, this is unstable and should switch to log1p formulation
  logres-= log(1 + pow( (x-mean)/shape ,2));
  if(give_log) return logres; else return exp(logres);
}


// Half normal distribution
   template<class Type>
   Type ldhalfnorm(Type x, Type var){
   return 0.5*log(2)-0.5*log(var*M_PI)+pow(x,2)/(2*var);
   }


// helper function to make sparse SPDE precision matrix
// Inputs:
// logkappa: log(kappa) parameter value
// logtau: log(tau) parameter value
//  M0, M1, M2: these sparse matrices are output from R::INLA::inla.spde2.matern()$param.inla$M*
template<class Type>
  SparseMatrix<Type> spde_Q(Type logkappa, Type logtau, SparseMatrix<Type> M0, SparseMatrix<Type> M1, SparseMatrix<Type> M2) {
  SparseMatrix<Type> Q;
  Type kappa2 = exp(2. * logkappa);
  Type kappa4 = kappa2*kappa2;
  Q = pow(exp(logtau), 2.)  * (kappa4*M0 + Type(2.0)*kappa2*M1 + M2);
  return Q;
}
template<class Type>
Type objective_function<Type>::operator() ()
{


//=============================================================================================================
//                                              DATA SECTION
//=============================================================================================================
  
// Vectors of real data
   DATA_VECTOR(y_i);         // y_i (response variable)
   //DATA_VECTOR(x_1);         // covariate
   

// SPDE objects
   DATA_SPARSE_MATRIX(M0);    // used to make gmrf precision
   DATA_SPARSE_MATRIX(M1);    // used to make gmrf precision
   DATA_SPARSE_MATRIX(M2);    // used to make gmrf precision
   
   DATA_SPARSE_MATRIX(Aproj); // used to project from spatial mesh to data locations
   
// Prior means
   DATA_SCALAR(tau0);  // mean of prior for logtauO
   DATA_SCALAR(kappa0); // mean of prior for logkappa
  
  
//=============================================================================================================
//                                              PARAMETERS SECTION
//=============================================================================================================
  
// Fixed effects
  PARAMETER(beta0);		                      // global mean
  //PARAMETER(beta1);		                      // beta 1
  PARAMETER(logsigma);		         
  PARAMETER(logtau);		                      // spatial process
  PARAMETER(logkappa);		                   // decorrelation distance (kind of)
  
  PARAMETER_VECTOR(omega_s);	               // spatial effects
  
  SparseMatrix<Type> Q   = spde_Q(logkappa, logtau, M0, M1, M2);
  
  

//===================================
//               Priors
//===================================
   Type nlp=0.0;                          // negative log prior  (priors)
   
// beta0
   nlp -= dnorm(beta0, Type(0.0),   Type(1.0), true); // beta0

// beta1 
   //nlp -= dnorm(beta1, Type(0.0),   Type(1.0), true); // beta0
   
// Variance component
   Type sigma = exp(logsigma);
   nlp -= dcauchy(sigma,   Type(0.0), Type(5.0));
  
  
// Hyperpriors
   Type tau    = exp(logtau);
   Type kappa  = exp(logkappa);
   
   nlp -= dnorm(tau,    tau0,     Type(1.0), true);
   nlp -= dnorm(kappa,  kappa0,   Type(1.0), true);
   
   
   
//=============================================================================================================
// Objective function is sum of negative log likelihood components
   using namespace density;
   int n_i = y_i.size();	             // number of observations 
   vector<Type> projS(n_i);             // value of gmrf at data points
   Type nll_omega=0;		               // spatial effects
  
   
// Probability of random effects
   nll_omega += SCALE(GMRF(Q), 1/exp(logtau) )(omega_s);
   //nll_omega += GMRF(Q)(omega_s);
   
// The model predicted for each observation, in natural space:
   vector<Type> mu(n_i);
   projS = Aproj * omega_s; // Project S at mesh pts to data pts
   for( int i=0; i<n_i; i++){
      //mu(i) = beta0 + beta1*x_1(i) + projS(i);
      mu(i) = beta0 + projS(i);
}
  
    
// Probability of the data, given random effects (likelihood)
   vector<Type> log_lik(n_i);
   for( int i = 0; i<n_i; i++){
       log_lik(i) = dnorm(y_i(i), mu(i), sigma, true);
}
  
   Type log_like = -log_lik.sum(); // total NLL
   

// Simule data from mu
   vector<Type> y_sim(n_i);
   for( int i=0; i<n_i; i++){
      SIMULATE {
      y_sim(i) = rnorm(mu(i), sigma);
    REPORT(y_sim)};
}


// Jacobian adjustment 
   nlp -= logsigma + logtau + logkappa;
   
   
// Calculate joint negative log likelihood
   Type jnll = log_like + nll_omega + nlp;
  
   vector<Type> preds = mu;
    
    
// Derived quantities, given parameters
// Geospatial
   Type rho = sqrt(8.0) / exp(logkappa);
   //Type sigma_s = 1 / sqrt(4 * M_PI * 2*tau * 2*kappa);
   //Type sigma_s = 1 / sqrt(4 * M_PI) * tau * kappa);
   Type sigma_s = 1.0 / sqrt (4.0 * M_PI*exp (2.0 * logtau ) * exp (2.0 * logkappa ));
   
   
//=====================================================================================================
// Reporting
   REPORT(jnll);
   REPORT(log_like);
   REPORT(nll_omega);
  
   REPORT(beta0);
   //REPORT(beta1);
   REPORT(omega_s);
   REPORT(preds);
   REPORT(logsigma);
   REPORT(logtau);
   REPORT(logkappa);
   REPORT(log_lik);
   REPORT(sigma);
   REPORT(tau);
   REPORT(kappa);
   REPORT(rho);		         
   REPORT(sigma_s);		         
  
//=======================================================================================================
// AD report (standard devations)
   ADREPORT(beta0);	                  
   //ADREPORT(beta1);	
   ADREPORT(logsigma);
   ADREPORT(logtau);
   ADREPORT(logkappa);
    
// Derived geospatial components
   ADREPORT(sigma);
   ADREPORT(rho);		               // geostatistical range
   ADREPORT(sigma_s);		         
   ADREPORT(tau);		               
   ADREPORT(kappa);		         
   return jnll;
}


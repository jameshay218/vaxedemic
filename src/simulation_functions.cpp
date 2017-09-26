#include <RcppArmadillo.h>
#include <RcppEigen.h>
using namespace std;
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

using namespace Eigen;

//' Matrix multiplication test
//' 
//' Testing
//' @param x
//' @param y
//' @return the multiplied matrix
//' @export
//' @useDynLib
//[[Rcpp::export]]
arma::mat matrix_mult_armadillo(arma::mat x, arma::mat y){
  return(x*y);
}


//' Matrix multiplication
//' 
//' Uses RcppEigen to multiply two matrices quickly
//' @param x the first matrix
//' @param y the second matrix
//' @return the multiplied matrix
//' @export
//' @useDynLib
//[[Rcpp::export]]
Eigen::MatrixXd matrix_mult(const Eigen::MatrixXd &x, const Eigen::MatrixXd &y){
  return(x*y);
}

double convert_poisson(double _lambda){
  return(1 - exp(-_lambda));
}

//' Simulation
//'
//' Wow
//' @export
//[[Rcpp::export]]
Rcpp::NumericMatrix stochastic_simulation(
					  Eigen::MatrixXd S,
					  Eigen::MatrixXd I,
					  Eigen::MatrixXd R,
					  double beta,
					  double gamma,
					  Eigen::MatrixXd Kdelta,
					  Eigen::MatrixXd K1,
					  Eigen::MatrixXd Mm1,
					  int div,
					  int ndays,
					  int n_age_groups,
					  int n_countries){

  MatrixXd Smat(n_age_groups*n_countries,ndays*div);
  MatrixXd Imat(n_age_groups*n_countries,ndays*div);
  MatrixXd tmp(I.rows(), I.cols());
  MatrixXd tmp1(I.rows(), I.cols());
  
  MatrixXd lambda(n_countries,1);
  MatrixXd Pinf(n_countries,1);

  double newInfections = 0;
  double newRecoveries = 0;
  
  int nsteps = ndays*div;
  
  for(int t = 0; t < ndays*div; ++t){
    tmp = K1*I;
    tmp1 = tmp.array() * Mm1.array();
    lambda = beta*Kdelta*tmp1;
    Pinf = lambda.unaryExpr(&convert_poisson);
    
    // Calculate new infections
    for(int i = 0; i < S.rows(); ++i){
      newInfections = R::rbinom(S(i,0),Pinf(i,0));
      I(i,0) += newInfections;
      S(i,0) -= newInfections;
    }

    // Calculate new recoveries
    for(int i = 0; i < I.rows(); ++i){
      newRecoveries = R::rbinom(I(i,0),1-exp(-gamma));
      I(i, 0) -= newRecoveries;
    }
    
    Imat.col(t) = I.col(0);

  }
  return(Rcpp::wrap(Imat));


  /*
 for(int i=0; i < nsteps; ++i){
lamb
    lambda = matrix_mult(beta*Kdelta, matrix_mult(K1, I)*Mm1);
    Rcpp::Rcout << lambda(0) << std::endl;

    for_each(lambda.begin(), lambda.end(), convert_poisson);
    Pinf = 1 - exp(-lambda);

    newInfections = R::rbinom(n_age_groups*n_countries, S, Pinf)
  */
  

  //return(Rcpp::wrap(Smat));
}
				    
				    

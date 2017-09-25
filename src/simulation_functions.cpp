#include <RcppArmadillo.h>
#include <RcppEigen.h>
using namespace std;

// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

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
Eigen::MatrixXd matrix_mult(Eigen::MatrixXd x, Eigen::MatrixXd y){
  return(x*y);
}

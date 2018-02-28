#include <Rcpp.h>
using namespace Rcpp;

//' Placeholder
//'
//[[Rcpp::export]]
double convert_poisson(double lambda){
  return(1 - exp(-lambda));
}

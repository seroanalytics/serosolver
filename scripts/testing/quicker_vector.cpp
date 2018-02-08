#include <Rcpp.h>
using namespace Rcpp;

#define MAX(a,b) ((a) < (b) ? (b) : (a)) // define MAX function for use later

// [[Rcpp::export]]
NumericVector f_vector(NumericVector x, double sigma) {
  NumericVector x2(x.size());
  for(int i = 0; i < x2.size(); ++i){
    x2[i] = MAX(1 - x[i]*sigma,0);
  }
  return(x2);
}

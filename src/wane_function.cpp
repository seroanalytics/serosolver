#include "wane_function.h"

//' Function to calculate non-linear waning
//'  All additional parameters for the function are declared here
//' @param theta NumericVector, the named vector of model parameters
//' @param time_infected double the time infected (samplingTime - circulation_time)
//' @return value of waning parameter based on time since infected 
//' @useDynLib serosolver
//' @export
// [[Rcpp::export]]
double wane_function(NumericVector theta, double time_infected, double wane){
  // Declare variables 
  double kappa = theta["kappa"];
  double t_change = theta["t_change"];
  double wane_2 = -kappa*wane;
  double wane_2_val; // Interaction term

  // Calculate the interaction term
  if(time_infected > t_change){
    wane_2_val = wane_2*(time_infected - t_change); 
  }else{
    wane_2_val = 0;
  }
  return (wane*time_infected+wane_2_val);
}

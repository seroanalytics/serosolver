#include "likelihood_funcs.h"

#include <iostream>
#include <fstream>
//' Marginal prior probability (p(Z)) of a particular infection history matrix single prior
//'  Prior is independent contribution from each year
//' @param infection_history IntegerMatrix, the infection history matrix
//' @param n_alive IntegerVector, vector giving the number of individuals alive in each time unit
//' @param shape1 double, shape1 (alpha) parameter for beta distribution prior
//' @param shape2 double, shape2 (beta) parameter for beta distribution prior
//' @return a single prior probability
//' @export
//' @family inf_mat_prior
// [[Rcpp::export]]
double inf_mat_prior_cpp(const IntegerMatrix& infection_history, const IntegerVector& n_alive, double shape1, double shape2){
  // Prior on each year
  double m, n;
  double lik=0;
  double lbeta_const = R::lbeta(shape1, shape2); // Contribution of prior
  for(int i = 0; i < n_alive.size(); ++i){
    m = sum(infection_history(_,i)); // Number of infections in that year
    n = n_alive(i); // Number of individuals alive in that year
    lik += R::lbeta(m+shape1,n-m+shape2)-lbeta_const; // Contribution of augmented data and prior for that year
  }
  return(lik);
}

//' Marginal prior probability (p(Z)) of a particular infection history matrix vector prior
//'  Prior is independent contribution from each year, but each year has its own shape parameters
//' @param infection_history IntegerMatrix, the infection history matrix
//' @param n_alive IntegerVector, vector giving the number of individuals alive in each time unit
//' @param shape1s NumericVector, shape1 (alpha) parameters for beta distribution prior, one for each time unit
//' @param shape2s NumericVector, shape2 (beta) parameters for beta distribution prior, one for each time unit
//' @return a single prior probability
//' @export
//' @family inf_mat_prior
// [[Rcpp::export]]
double inf_mat_prior_cpp_vector(const IntegerMatrix& infection_history, const IntegerVector& n_alive, 
                                const NumericVector& shape1s, const NumericVector& shape2s){
  // Prior on each year
  double m, n;
  double lik=0;
  for(int i = 0; i < n_alive.size(); ++i){ 
    m = sum(infection_history(_,i)); // Number of infections in that year
    n = n_alive(i); // Number of individuals alive in that year
    if(n > 0){
      lik += R::lbeta(m+shape1s[i],n-m+shape2s[i]) - 
        R::lbeta(shape1s[i], shape2s[i]); // Contribution of augmented data and prior for that year
    }
  }
  return(lik);
}

//' Marginal prior probability (p(Z)) of infection history matrix with groups
//'  Prior is independent contribution from each year and group
//' @param n_infections IntegerMatrix, the total number of infections in each time point/group
//' @param n_alive IntegerMatrix, matrix giving the number of individuals alive in each time unit in each group
//' @param shape1 double, shape1 (alpha) parameter for beta distribution prior
//' @param shape2 double, shape2 (beta) parameter for beta distribution prior
//' @return a single prior probability
//' @export
//' @family inf_mat_prior
// [[Rcpp::export]]
double inf_mat_prior_group_cpp(const IntegerMatrix& n_infections, const IntegerMatrix& n_alive, double shape1, double shape2){
  // Prior on each time
  double m, n;
  double lik=0;
  double lbeta_const = R::lbeta(shape1, shape2); // Contribution of prior
  for(int j = 0; j < n_alive.nrow(); ++j){
    for(int i = 0; i < n_alive.ncol(); ++i){
      m = n_infections(j,i);
      n = n_alive(j,i); // Number of individuals alive in that time
      // Only solve if n > 0
      if(n > 0){
	lik += R::lbeta(m+shape1,n-m+shape2)-lbeta_const; // Contribution of augmented data and prior for that time
      }
    }
  }
  return(lik);
}

//' Marginal prior probability (p(Z)) of a particular infection history matrix vector prior by group
//'  Prior is independent contribution from each time, but each time has its own shape parameters
//' @param n_infections IntegerMatrix, the total number of infections in each time point/group
//' @param n_alive IntegerMatrix, matrix giving the number of individuals alive in each time unit in each group
//' @param shape1s NumericVector, shape1 (alpha) parameters for beta distribution prior, one for each time unit
//' @param shape2s NumericVector, shape2 (beta) parameters for beta distribution prior, one for each time unit
//' @return a single prior probability
//' @export
//' @family inf_mat_prior
// [[Rcpp::export]]
double inf_mat_prior_group_cpp_vector(const IntegerMatrix& n_infections, const IntegerMatrix& n_alive, 
                                      const NumericVector& shape1s, const NumericVector& shape2s){
  // Prior on each time
  double m, n;
  double lik=0;
  for(int j = 0; j < n_alive.nrow(); ++j){
    for(int i = 0; i < n_alive.ncol(); ++i){ 
      m = n_infections(j,i); // Number of infections in that time
      n = n_alive(j,i); // Number of individuals alive in that time
      lik += R::lbeta(m+shape1s[i],n-m+shape2s[i]) - R::lbeta(shape1s[i], shape2s[i]); // Contribution of augmented data and prior for that time
    }
  }
  return(lik);
}

//' Marginal prior probability (p(Z)) of a particular infection history matrix total prior
//'  Prior here is on the total number of infections across all individuals and times
//' @param n_infections_group IntegerVector, the total number of infections in each group
//' @param n_alive_group IntegerVector, vector giving total number of potential infection events per group
//' @param shape1 double, shape1 (alpha) parameter for beta distribution prior
//' @param shape2 double, shape2 (beta) parameter for beta distribution prior
//' @return a single prior probability
//' @export
//' @family inf_mat_prior
// [[Rcpp::export]]
double inf_mat_prior_total_group_cpp(const IntegerVector& n_infections_group, const IntegerVector& n_alive_group, double shape1, double shape2){
  double lik=0;
  double beta_const = R::lbeta(shape1, shape2);
  int n_infections =0;
  for(int i = 0; i < n_alive_group.size(); ++i){
    n_infections = n_infections_group(i);
    lik += R::lbeta(n_infections + shape1, n_alive_group(i) - n_infections + shape2) - beta_const;
  }
  return(lik);
}


//' Fast observation error function
//'  Calculate the probability of a set of observed antibody levels given a corresponding set of predicted antibody levels. 
//' @param theta NumericVector, a named parameter vector giving the normal distribution standard deviation and the max observable antibody level
//' @param obs NumericVector, the vector of observed log antibody levels
//' @param predicted_antibody_levels NumericVector, the vector of predicted log antibody levels
//' @param a vector of same length as the input data giving the probability of observing each observation given the predictions
//' @return a likelihood for each observed antibody level
//' @export
//' @family likelihood_functions
// [[Rcpp::export(rng = false)]]
NumericVector likelihood_func_fast(const NumericVector &theta, const NumericVector &obs, const NumericVector &predicted_antibody_levels){
  int total_measurements = predicted_antibody_levels.size();
  NumericVector ret(total_measurements);
  const double sd = theta["obs_sd"];
  const double den = sd*M_SQRT2;
  const double max_measurement = theta["max_measurement"];
  const double min_measurement = theta["min_measurement"];
  const double log_const = log(0.5);

  for(int i = 0; i < total_measurements; ++i){
    // Most antibody levels are between 0 and max_measurement, this is the difference in normal cdfs
    if(obs[i] < max_measurement && obs[i] >= (min_measurement + 1)){
      ret[i] = log_const + log((erf((obs[i] + 1.0 - predicted_antibody_levels[i]) / den) -
				erf((obs[i]     - predicted_antibody_levels[i]) / den)));    
      // For antibody levels above the maximum, 
    } else if(obs[i] >= max_measurement) {
      ret[i] = log_const + log(erfc((max_measurement - predicted_antibody_levels[i])/den));
    } else {
      ret[i] = log_const + log(1.0 + erf(((min_measurement + 1) - predicted_antibody_levels[i])/den));
    }
  }
  return(ret);
} 


//' Fast observation error function continuous
//'  Calculate the probability of a set of observed antibody levels given a corresponding set of predicted antibody levels assuming continuous, bounded observations.
//' @name Fast observation error function continuous
//' @param theta NumericVector, a named parameter vector giving the normal distribution standard deviation and the max observable antibody level
//' @param obs NumericVector, the vector of observed log antibody levels
//' @param predicted_antibody_levels NumericVector, the vector of predicted log antibody levels
//' @param a vector of same length as the input data giving the probability of observing each observation given the predictions
//' @return a likelihood for each observed antibody level
//' @export
//' @family likelihood_functions
// [[Rcpp::export(rng = false)]]
NumericVector likelihood_func_fast_continuous(const NumericVector &theta, const NumericVector &obs, const NumericVector &predicted_antibody_levels){
 int total_measurements = predicted_antibody_levels.size();
 NumericVector ret(total_measurements);
 const double sd = theta["obs_sd"];
 const double den = sd*M_SQRT2;
 const double den2 = log(sd*2.50662827463);
 const double max_measurement = theta["max_measurement"];
 const double min_measurement = theta["min_measurement"];
 const double log_const = log(0.5);

 for(int i = 0; i < total_measurements; ++i){
   // Most antibody levels are between 0 and max_measurement, this is the difference in normal cdfs
   if(obs[i] < max_measurement && obs[i] > min_measurement){
     ret[i] = -0.5*(pow((obs[i]-predicted_antibody_levels[i])/sd, 2)) - den2;
     // For antibody levels above the maximum, 
   } else if(obs[i] >= max_measurement) {
     ret[i] = log_const + log(erfc((max_measurement - predicted_antibody_levels[i])/den));
   } else {
     ret[i] = log_const + log(1.0 + erf((min_measurement - predicted_antibody_levels[i])/den));
   }
 }
 return(ret);
} 

// Likelihood calculation for infection history proposal
// Not really to be used elsewhere other than in \code{\link{inf_hist_prop_prior_v2_and_v4}}, as requires correct indexing for the predicted antibody levels vector. Also, be very careful, as predicted_antibody_levels is set to 0 at the end!
void proposal_likelihood_func(double &new_prob,
			      NumericVector &predicted_antibody_levels,
			      const int &indiv,
			      const NumericVector &data,
			      const NumericVector &repeat_data,
			      const IntegerVector &repeat_indices,
			      const IntegerVector &cum_nrows_per_individual_in_data,
			      const IntegerVector &cum_nrows_per_individual_in_repeat_data,
			      const double &log_const,
			      const double &den,
			      const double &max_measurement,
			      const double &min_measurement,
			      const bool &repeat_data_exist,
			      const double &obs_weight = 1.0){
  
  for(int x = cum_nrows_per_individual_in_data[indiv]; x < cum_nrows_per_individual_in_data[indiv+1]; ++x){
    if(data[x] < max_measurement && data[x] >= (min_measurement + 1)){
      new_prob += log_const + log((erf((data[x] + 1.0 - predicted_antibody_levels[x]) / den) -
				   erf((data[x]     - predicted_antibody_levels[x]) / den)));    
    } else if(data[x] >= max_measurement) {
      new_prob += log_const + log(erfc((max_measurement - predicted_antibody_levels[x])/den));
    } else {
      new_prob += log_const + log(1.0 + erf(((min_measurement + 1) - predicted_antibody_levels[x])/den));
    }
  }

  // =====================
  // Do something for repeat data here
  if(repeat_data_exist){
    for(int x = cum_nrows_per_individual_in_repeat_data[indiv]; x < cum_nrows_per_individual_in_repeat_data[indiv+1]; ++x){
      if(repeat_data[x] < max_measurement && repeat_data[x] >= (min_measurement + 1)){
	new_prob += log_const + log((erf((repeat_data[x] + 1.0 - predicted_antibody_levels[repeat_indices[x]]) / den) -
				     erf((repeat_data[x]     - predicted_antibody_levels[repeat_indices[x]]) / den)));    
      } else if(repeat_data[x] >= max_measurement) {
	new_prob += log_const + log(erfc((max_measurement - predicted_antibody_levels[repeat_indices[x]])/den));
      } else {
	new_prob += log_const + log(1.0 + erf(((min_measurement + 1) - predicted_antibody_levels[repeat_indices[x]])/den));
      }
    }
  }
  
  // Re-weight likelihood
  new_prob = new_prob*obs_weight;
  
  // Need to erase the predicted antibody level data...
  for(int x = cum_nrows_per_individual_in_data[indiv]; x < cum_nrows_per_individual_in_data[indiv+1]; ++x){
    predicted_antibody_levels[x] = 0;
  }
}



// Likelihood calculation for infection history proposal for continuous, bounded data
// Not really to be used elsewhere other than in \code{\link{inf_hist_prop_prior_v2_and_v4}}, as requires correct indexing for the predicted antibody levels vector. Also, be very careful, as predicted_antibody_levels is set to 0 at the end!
void proposal_likelihood_func_continuous(double &new_prob,
                              NumericVector &predicted_antibody_levels,
                              const int &indiv,
                              const NumericVector &data,
                              const NumericVector &repeat_data,
                              const IntegerVector &repeat_indices,
                              const IntegerVector &cum_nrows_per_individual_in_data,
                              const IntegerVector &cum_nrows_per_individual_in_repeat_data,
                              const double &log_const,
                              const double &sd,
                              const double &den,
                              const double &den2,
                              const double &max_measurement,
                              const double &min_measurement,
                              const bool &repeat_data_exist,
                              const double &obs_weight = 1.0){
  for(int x = cum_nrows_per_individual_in_data[indiv]; x < cum_nrows_per_individual_in_data[indiv+1]; ++x){
    if(data[x] < max_measurement && data[x] > min_measurement){
      new_prob += -0.5*(pow((data[x]-predicted_antibody_levels[x])/sd, 2)) - den2;
    } else if(data[x] >= max_measurement) {
      new_prob += log_const + log(erfc((max_measurement - predicted_antibody_levels[x])/den));
    } else {
      new_prob += log_const + log(1.0 + erf((min_measurement - predicted_antibody_levels[x])/den));
    }
  }
  
  // =====================
  // Do something for repeat data here
  if(repeat_data_exist){
    for(int x = cum_nrows_per_individual_in_repeat_data[indiv]; x < cum_nrows_per_individual_in_repeat_data[indiv+1]; ++x){
      if(repeat_data[x] < max_measurement && repeat_data[x] > min_measurement){
        new_prob += -0.5*(pow((repeat_data[x]-predicted_antibody_levels[repeat_indices[x]])/sd, 2)) - den2;
      } else if(repeat_data[x] >= max_measurement) {
        new_prob += log_const + log(erfc((max_measurement - predicted_antibody_levels[repeat_indices[x]])/den));
      } else {
        new_prob += log_const + log(1.0 + erf((min_measurement - predicted_antibody_levels[repeat_indices[x]])/den));
      }
    }
  }
  // Re-weight likelihood
  new_prob = new_prob*obs_weight;
  
  // Need to erase the predicted antibody level data...
  for(int x = cum_nrows_per_individual_in_data[indiv]; x < cum_nrows_per_individual_in_data[indiv+1]; ++x){
    predicted_antibody_levels[x] = 0;
  }
}

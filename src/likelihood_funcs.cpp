#include "likelihood_funcs.h"

//' Marginal prior probability (p(Z)) of a particular infection history matrix single prior
//'  Prior is independent contribution from each year
//' @param infection_history IntegerMatrix, the infection history matrix
//' @param n_alive IntegerVector, vector giving the number of individuals alive in each time unit
//' @param alpha double, alpha parameter for beta distribution prior
//' @param beta double, beta parameter for beta distribution prior
//' @return a single prior probability
//' @export
//' @family inf_mat_prior
// [[Rcpp::export]]
double inf_mat_prior_cpp(const IntegerMatrix& infection_history, const IntegerVector& n_alive, double alpha, double beta){
  // Prior on each year
  double m, n;
  double lik=0;
  double lbeta_const = R::lbeta(alpha, beta); // Contribution of prior
  for(int i = 0; i < n_alive.size(); ++i){
    m = sum(infection_history(_,i)); // Number of infections in that year
    n = n_alive(i); // Number of individuals alive in that year
    lik += R::lbeta(m+alpha,n-m+beta)-lbeta_const; // Contribution of augmented data and prior for that year
  }
  return(lik);
}

//' Marginal prior probability (p(Z)) of a particular infection history matrix vector prior
//'  Prior is independent contribution from each year, but each year has its own alpha and beta
//' @param infection_history IntegerMatrix, the infection history matrix
//' @param n_alive IntegerVector, vector giving the number of individuals alive in each time unit
//' @param alphas NumericVector, alpha parameters for beta distribution prior, one for each time unit
//' @param betas NumericVector, beta parameters for beta distribution prior, one for each time unit
//' @return a single prior probability
//' @export
//' @family inf_mat_prior
// [[Rcpp::export]]
double inf_mat_prior_cpp_vector(const IntegerMatrix& infection_history, const IntegerVector& n_alive, const NumericVector& alphas, const NumericVector& betas){
  // Prior on each year
  double m, n;
  double lik=0;
  for(int i = 0; i < n_alive.size(); ++i){ 
    m = sum(infection_history(_,i)); // Number of infections in that year
    n = n_alive(i); // Number of individuals alive in that year
    lik += R::lbeta(m+alphas[i],n-m+betas[i]) - R::lbeta(alphas[i], betas[i]); // Contribution of augmented data and prior for that year
  }
  return(lik);
}



//' Marginal prior probability (p(Z)) of a particular infection history matrix total prior
//'  Prior here is on the total number of infections across all individuals and years
//' @param infection_history IntegerMatrix, the infection history matrix
//' @param n_alive IntegerVector, vector giving the number of individuals alive in each year
//' @param alpha double, alpha parameter for beta distribution prior
//' @param beta double, beta parameter for beta distribution prior
//' @return a single prior probability
//' @export
//' @family inf_mat_prior
// [[Rcpp::export]]
double inf_mat_prior_total_cpp(const IntegerMatrix& infection_history, const int& n_alive, double alpha, double beta){
  double m, n;
  double lik=0;
  int n_infections = sum(infection_history);
  lik = R::lbeta(n_infections + alpha, n_alive - n_infections + beta) - R::lbeta(alpha, beta);
  return(lik);
}


//' Fast observation error function
//'  Calculate the probability of a set of observed titres given a corresponding set of predicted titres. FAST IMPLEMENTATION
//' @param theta NumericVector, a named parameter vector giving the normal distribution standard deviation and the max observable titre
//' @param obs NumericVector, the vector of observed log titres
//' @param predicted_titres NumericVector, the vector of predicted log titres
//' @param a vector of same length as the input data giving the probability of observing each observation given the predictions
//' @return a likelihood for each observed titre
//' @export
//' @family likelihood_functions
// [[Rcpp::export(rng = false)]]
NumericVector likelihood_func_fast(const NumericVector &theta, const NumericVector &obs, const NumericVector &predicted_titres){
  int total_titres = predicted_titres.size();
  NumericVector ret(total_titres);
  const double sd = theta["error"];
  const double den = sd*M_SQRT2;
  const double max_titre = theta["MAX_TITRE"];
  const double log_const = log(0.5);

  for(int i = 0; i < total_titres; ++i){
    // Most titres are between 0 and max_titre, this is the difference in normal cdfs
    if(obs[i] < max_titre && obs[i] >= 1.0){
      ret[i] = log_const + log((erf((obs[i] + 1.0 - predicted_titres[i]) / den) -
				erf((obs[i]     - predicted_titres[i]) / den)));    
      // For titres above the maximum, 
    } else if(obs[i] >= max_titre) {
      ret[i] = log_const + log(erfc((max_titre - predicted_titres[i])/den));
    } else {
      ret[i] = log_const + log(1.0 + erf((1.0 - predicted_titres[i])/den));
    }
  }
  return(ret);
} 

// Likelihood calculation for infection history proposal
// Not really to be used elsewhere other than in \code{\link{infection_history_proposal_gibbs_fast}}, as requires correct indexing for the predicted titres vector. Also, be very careful, as predicted_titres is set to 0 at the end!
void proposal_likelihood_func(double &new_prob,
			      NumericVector &predicted_titres,
			      const int &indiv,
			      const NumericVector &data,
			      const NumericVector &repeat_data,
			      const IntegerVector &repeat_indices,
			      const IntegerVector &cum_nrows_per_individual_in_data,
			      const IntegerVector &cum_nrows_per_individual_in_repeat_data,
			      const double &log_const,
			      const double &den,
			      const double &max_titre){
  for(int x = cum_nrows_per_individual_in_data[indiv]; x < cum_nrows_per_individual_in_data[indiv+1]; ++x){
    if(data[x] < max_titre && data[x] >= 1.0){
      new_prob += log_const + log((erf((data[x] + 1.0 - predicted_titres[x]) / den) -
				   erf((data[x]     - predicted_titres[x]) / den)));    
    } else if(data[x] >= max_titre) {
      new_prob += log_const + log(erfc((max_titre - predicted_titres[x])/den));
    } else {
      new_prob += log_const + log(1.0 + erf((1.0 - predicted_titres[x])/den));
    }
  }

  // =====================
  // Do something for repeat data here
  for(int x = cum_nrows_per_individual_in_repeat_data[indiv]; x < cum_nrows_per_individual_in_repeat_data[indiv+1]; ++x){
    if(repeat_data[x] < max_titre && repeat_data[x] >= 1.0){
      new_prob += log_const + log((erf((repeat_data[x] + 1.0 - predicted_titres[repeat_indices[x]]) / den) -
				   erf((repeat_data[x]     - predicted_titres[repeat_indices[x]]) / den)));    
    } else if(repeat_data[x] >= max_titre) {
      new_prob += log_const + log(erfc((max_titre - predicted_titres[repeat_indices[x]])/den));
    } else {
      new_prob += log_const + log(1.0 + erf((1.0 - predicted_titres[repeat_indices[x]])/den));
    }
  }
  // Need to erase the predicted titre data...
  for(int x = cum_nrows_per_individual_in_data[indiv]; x < cum_nrows_per_individual_in_data[indiv+1]; ++x){
    predicted_titres[x] = 0;
  }
}




//' Calculate likelihood basic
//'
//' Calculates the likelihood of a given set of observed titres given predicted titres. Based on truncated discritised normal. DEPRECATED, as this is a very slow (but obvious) implementation
//' @param expected NumericVector, as returned by \code{\link{infection_model_indiv}}
//' @param data NumericVector, the vector of observed titres
//' @param theta NumericVector, the vector of named model parameters, requiring MAX_TITRE and error
//' @param titre_shifts NumericVector, OPTIONAL if using measurement bias, gives the shift to add to each expected titre
//' @return a single log likelihood value
//' @export
//[[Rcpp::export]]
double likelihood_titre_basic(const NumericVector &expected, 
			const NumericVector &data, 
			const NumericVector &theta,
			const NumericVector &titre_shifts
			){

  NumericVector expected1 = expected;
  if(titre_shifts.size() == expected.size()){
    expected1 = expected1 + titre_shifts;
  }
  int n = expected.size();
  double lnlike = 0;
  
  for(int i=0; i < n; ++i){
    if(!NumericVector::is_na(data[i])){
      if(data[i] >= theta["MAX_TITRE"]){
	lnlike += R::pnorm(theta["MAX_TITRE"], expected[i], theta["error"],0,1);
      } else if(data[i] < 1.0){
	lnlike += R::pnorm(1, expected[i], theta["error"],1,1);
      } else {
	lnlike += log(R::pnorm(data[i]+1, expected[i],theta["error"], 1,0) - 
		      R::pnorm(data[i], expected[i], theta["error"],1,0));
      }
    }
  }
  return(lnlike);
}

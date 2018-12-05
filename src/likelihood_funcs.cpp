#include "likelihood_funcs.h"

//' Marginal prior probability (p(Z)) of a particular infection history matrix
//' 
//' @param infHist IntegerMatrix, the infection history matrix
//' @param n_alive IntegerVector, vector giving the number of individuals alive in each time unit
//' @param alpha NumericVector, alpha parameter for beta distribution prior, one for each time unit
//' @param beta NumericVector, beta parameter for beta distribution prior, one for each time unit
//' @export
// [[Rcpp::export]]
double inf_mat_prior_cpp(const IntegerMatrix& infHist, const IntegerVector& n_alive, double alpha, double beta){
//double inf_mat_prior_cpp(const IntegerMatrix& infHist, const IntegerVector& n_alive, const NumericVector& alphas, const NumericVector& betas){
  // Prior on each year
  double m, n;
  double lik=0;
  double lbeta_const = R::lbeta(alpha, beta);
  for(int i = 0; i < n_alive.size(); ++i){
    m = sum(infHist(_,i));
    n = n_alive(i);
    lik += R::lbeta(m+alpha,n-m+beta)-lbeta_const;
  }
  return(lik);
}


//' Marginal prior probability (p(Z)) of a particular infection history matrix
//' 
//' @param infHist IntegerMatrix, the infection history matrix
//' @param n_alive IntegerVector, vector giving the number of individuals alive in each year
//' @param alpha double, alpha parameter for beta distribution prior
//' @param beta double, beta parameter for beta distribution prior
//' @export
// [[Rcpp::export]]
double inf_mat_prior_total_cpp(const IntegerMatrix& infHist, const int& n_alive, double alpha, double beta){
  double m, n;
  double lik=0;
  int n_infections = sum(infHist);
  lik = R::lbeta(n_infections + alpha, n_alive - n_infections + beta) - R::lbeta(alpha, beta);
  return(lik);
}


//' @export
// [[Rcpp::export(rng = false)]]
NumericVector likelihood_func_fast(NumericVector theta, NumericVector obs, NumericVector predicted_titres){
  int total_titres = predicted_titres.size();
  NumericVector ret(total_titres);
  const double sd = theta["error"];
  const double den = sd*M_SQRT2;
  const double max_titre = theta["MAX_TITRE"];
  const double log_const = log(0.5);

  for(int i = 0; i < total_titres; ++i){
    if(obs[i] <= max_titre && obs[i] > 0.0){
      ret[i] = log_const + log((erf((obs[i] + 1.0 - predicted_titres[i]) / den) -
				erf((obs[i]     - predicted_titres[i]) / den)));    
    } else if(obs[i] > max_titre) {
      ret[i] = log_const + log(1.0 + erf((max_titre - predicted_titres[i])/den));
    } else {
      ret[i] = log_const + log(1.0 + erf((1.0 - predicted_titres[i])/den));
    }
  }
  return(ret);
} 


double erf_native(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return sign*y;
}

//' @export
// [[Rcpp::export]]
NumericVector likelihood_func_fast_native(const NumericVector &theta, const NumericVector &obs, const NumericVector &predicted_titres){
  int total_titres = predicted_titres.size();
  NumericVector ret(total_titres);
  const double sd = theta["error"];
  const double den = sd*M_SQRT2;
  const double max_titre = theta["MAX_TITRE"];
  const double log_const = log(0.5);

  for(int i = 0; i < total_titres; ++i){
    if(obs[i] <= max_titre && obs[i] > 0.0){
      ret[i] = log_const + log((erf_native((obs[i] + 1.0 - predicted_titres[i]) / den) -
				erf_native((obs[i]     - predicted_titres[i]) / den)));    
    } else if(obs[i] > max_titre) {
      ret[i] = log_const + log(1.0 + erf_native((max_titre - predicted_titres[i])/den));
    } else {
      ret[i] = log_const + log(1.0 + erf_native((1.0 - predicted_titres[i])/den));
    }
  }
  return(ret);
} 


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
    if(data[x] <= max_titre && data[x] > 0.0){
      new_prob += log_const + log((erf((data[x] + 1.0 - predicted_titres[x]) / den) -
				   erf((data[x]     - predicted_titres[x]) / den)));    
    } else if(data[x] > max_titre) {
      new_prob += log_const + log(1.0 + erf((max_titre - predicted_titres[x])/den));
    } else {
      new_prob += log_const + log(1.0 + erf((1.0 - predicted_titres[x])/den));
    }
  }

  // =====================
  // Do something for repeat data here
  for(int x = cum_nrows_per_individual_in_repeat_data[indiv]; x < cum_nrows_per_individual_in_repeat_data[indiv+1]; ++x){
    if(repeat_data[x] <= max_titre && repeat_data[x] > 0.0){
      new_prob += log_const + log((erf((repeat_data[x] + 1.0 - predicted_titres[repeat_indices[x]]) / den) -
				   erf((repeat_data[x]     - predicted_titres[repeat_indices[x]]) / den)));    
    } else if(repeat_data[x] > max_titre) {
      new_prob += log_const + log(1.0 + erf((max_titre - predicted_titres[repeat_indices[x]])/den));
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
//' @return a single log likelihood
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
      if(data[i] > theta["MAX_TITRE"]){
	lnlike += R::pnorm(theta["MAX_TITRE"], expected[i], theta["error"],0,1);
      } else if(data[i] <= 0){
	lnlike += R::pnorm(1, expected[i], theta["error"],1,1);
      } else {
	lnlike += log(R::pnorm(data[i]+1, expected[i],theta["error"], 1,0) - 
		      R::pnorm(data[i], expected[i], theta["error"],1,0));
      }
    }
  }
  return(lnlike);
}

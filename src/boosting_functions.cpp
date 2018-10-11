#include <Rcpp.h>
#include "boosting_functions.h"
using namespace Rcpp;

#define MAX(a,b) ((a) < (b) ? (b) : (a)) // define MAX function for use later

// [[Rcpp::export]]
void multiple_infection_base_boosting(NumericVector &predicted_titres,
				    const NumericVector &theta,
				      const IntegerVector &cumu_infection_history,
				      const IntegerVector &masked_infection_history,
				      const IntegerVector &infection_map_indices, 
				    const IntegerVector &measurement_map_indices,
				    const NumericVector &antigenic_map_long, 
				    const NumericVector &antigenic_map_short, 
				    const NumericVector &waning,
				    const int &number_strains
				    ){
  int n_samples = measurement_map_indices.size();
  int max_infections = masked_infection_history.size();
  double mu = theta["mu"];
  double mu_short = theta["mu_short"];
  double tau = theta["tau"];

  for(int i = 0; i < max_infections; ++i){
    for(int k = 0; k < n_samples; ++k){
      if(masked_infection_history[i] > 0){
	predicted_titres[k] += MAX(0, 1.0 - tau * (cumu_infection_history[i]-1.0))*
	  ((mu * antigenic_map_long[measurement_map_indices[k] * number_strains + infection_map_indices[i]]) + 
	   (mu_short * antigenic_map_short[measurement_map_indices[k] * number_strains + infection_map_indices[i]]) * 
	   waning[i]);    
      }
    }
  }
}

// [[Rcpp::export]]
void multiple_infection_titre_dependent_boost(NumericVector &predicted_titres, 
					      NumericVector &monitored_titres,
					      const NumericVector &theta,
					      const NumericVector &infection_times,
					      const IntegerVector &cumu_infection_history,
					      const IntegerVector &masked_infection_history,
					      const IntegerVector &infection_map_indices,
					      const IntegerVector &measurement_map_indices,
					      const NumericVector &antigenic_map_long, 
					      const NumericVector &antigenic_map_short, 
					      const NumericVector &waning,
					      const int &number_strains
					      ){
  double circulation_time;
  double monitored_titre = 0;
  double long_boost = 0;
  double short_boost = 0;
  double boost = 0;

  double mu = theta["mu"];
  double mu_short = theta["mu_short"];
  double tau = theta["tau"];
  double gradient = theta["gradient"];
  double boost_limit = theta["boost_limit"];
  double wane = theta["wane"];
  
  int max_infections = infection_times.size();
  int n_samples = measurement_map_indices.size();

  for(int i = 0; i < max_infections; ++i){
    circulation_time = infection_times[i];

    if(masked_infection_history[i] > 0){
      for(int ii = i - 1; ii >= 0; --ii){
	if(masked_infection_history[ii] > 0){
	  long_boost = MAX(0, 1.0 - tau*(cumu_infection_history[ii] - 1.0)) * // Antigenic seniority
	    (mu * antigenic_map_long[infection_map_indices[i] * 
				     number_strains + infection_map_indices[ii]]);      
      
	  // Short term cross reactive boost
	  short_boost =  MAX(0, 1.0 - tau*(cumu_infection_history[ii] - 1.0)) * // Antigenic seniority
	    (mu_short * antigenic_map_short[infection_map_indices[i] * 
					    number_strains + infection_map_indices[ii]]);

	  if(monitored_titres[ii] >= boost_limit){
	    long_boost =  long_boost * (1 - gradient * boost_limit); // Titre dependent boosting - at ceiling
	    short_boost =  short_boost * (1 - gradient * boost_limit); // Titre dependent boosting - at ceiling
	  } else {
	    long_boost = long_boost * (1 - gradient * monitored_titres[ii]); // Titre dependent boosting - below ceiling
	    short_boost = short_boost * (1 - gradient * monitored_titres[ii]); // Titre dependent boosting - below ceiling
	  }
	  long_boost = MAX(0, long_boost);
	  short_boost = MAX(0, short_boost);
	  boost = long_boost + short_boost * MAX(0, 1.0 - wane * (circulation_time - infection_times[ii]));
	  monitored_titre += boost;
	}
	monitored_titres[i] = monitored_titre;
      }
      for(int k = 0; k < n_samples; ++k){
	// How much boosting experienced from this infection?
	long_boost = MAX(0, 1.0 - tau * (cumu_infection_history[i] - 1.0)) *
	  (mu * antigenic_map_long[measurement_map_indices[k] * 
				   number_strains + infection_map_indices[i]]);
    
	// Short term cross reactive boost
	short_boost = MAX(0, 1.0 - tau * (cumu_infection_history[i] - 1.0)) *
	  (mu_short * antigenic_map_short[measurement_map_indices[k] * 
					  number_strains + infection_map_indices[i]]);
    
    
	if(monitored_titres[i] >= boost_limit){
	  long_boost =  long_boost * (1 - gradient * boost_limit); // Titre dependent boosting - at ceiling
	  short_boost =  short_boost * (1 - gradient * boost_limit); // Titre dependent boosting - at ceiling
	} else {
	  long_boost = long_boost * (1 - gradient * monitored_titres[i]); // Titre dependent boosting - below ceiling
	  short_boost = short_boost * (1 - gradient * monitored_titres[i]); // Titre dependent boosting - below ceiling
	}
	long_boost = MAX(0, long_boost);
	short_boost = MAX(0, short_boost);
	boost = long_boost + short_boost * waning[i];
	predicted_titres[k] += boost;
      }
    }
  }
}

// [[Rcpp::export]]
void add_multiple_infections_boost(NumericVector &predicted_titres, 
				   NumericVector &monitored_titres,
				   const NumericVector &theta,
				   const NumericVector &infection_times,
				   const IntegerVector &cumu_infection_history,
				   const IntegerVector &masked_infection_history,
				   const IntegerVector &infection_map_indices,
				   const IntegerVector &measurement_map_indices,
				   const NumericVector &antigenic_map_long, 
				   const NumericVector &antigenic_map_short, 
				   const NumericVector &waning,
				   const int &number_strains,
				   const int &titre_dependent_boosting,
				   const int &age,
				   const Nullable<List> &additional_arguments){
  if(titre_dependent_boosting == 1){
    multiple_infection_titre_dependent_boost(predicted_titres, 
					     monitored_titres,
					     theta,
					     infection_times,
					     cumu_infection_history,
					     masked_infection_history,
					     infection_map_indices,
					     measurement_map_indices,
					     antigenic_map_long, 
					     antigenic_map_short, 
					     waning,
					     number_strains);    
  } else {
    multiple_infection_base_boosting(predicted_titres,
				     theta, 
				     cumu_infection_history,
				     masked_infection_history,
				     infection_map_indices,
				     measurement_map_indices,
				     antigenic_map_long,
				     antigenic_map_short,
				     waning,
				     number_strains);							       
  }
}
				

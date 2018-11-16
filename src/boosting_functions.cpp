#include <Rcpp.h>
#include "boosting_functions.h"
using namespace Rcpp;


#define MAX(a,b) ((a) < (b) ? (b) : (a)) // define MAX function for use later

//[[Rcpp::export]]
void multiple_infection_strain_dependent(NumericVector &predicted_titres,
					 const NumericVector &theta,
					 const IntegerVector &cumu_infection_history,
					 const IntegerVector &infection_map_indices, 
					 const IntegerVector &measurement_map_indices,
					 const NumericVector &antigenic_map_long, 
					 const NumericVector &antigenic_map_short, 
					 const NumericVector &waning,
					 const int &number_strains,
					 List additional_arguments){
  int n_samples = measurement_map_indices.size();
  int max_infections = cumu_infection_history.size();
  double mu;
  double mu_short = theta["mu_short"];
  double tau = theta["tau"];
  double n_inf, inf_map_index, wane;
  NumericVector mus = additional_arguments["mus"];
  IntegerVector boosting_vec_indices = additional_arguments["boosting_vec_indices"];
  for(int i = 0; i < max_infections; ++i){
    n_inf = cumu_infection_history[i] - 1.0;
    inf_map_index = infection_map_indices[i];
    mu = mus[boosting_vec_indices[inf_map_index]];
    //Rcpp::Rcout << "Mus: " << mus << std::endl;
    //Rcpp::Rcout << "Boosting vec index: " << boosting_vec_indices << std::endl;
    //Rcpp::Rcout << "Infection map index: " << inf_map_index << std::endl;
    //Rcpp::Rcout << "Mu used: " << mu << std::endl;
    wane = waning[i];
    //if(masked_infection_history[i] > 0){
    for(int k = 0; k < n_samples; ++k){
      predicted_titres[k] += MAX(0, 1.0 - tau * n_inf)*
	((mu * antigenic_map_long[measurement_map_indices[k] * number_strains + inf_map_index]) + 
	 (mu_short * antigenic_map_short[measurement_map_indices[k] * number_strains + inf_map_index]) * 
	 wane);    
      //  }
    }
  }  
}

//[[Rcpp::export]]
void multiple_infection_base_boosting(NumericVector &predicted_titres,
				      const NumericVector &theta,
				      const IntegerVector &cumu_infection_history,
				      const IntegerVector &infection_map_indices, 
				      const IntegerVector &measurement_map_indices,
				      const NumericVector &antigenic_map_long, 
				      const NumericVector &antigenic_map_short, 
				      const NumericVector &waning,
				      const NumericVector &seniority,
				      const int &number_strains,
				      const int &n_samples,
				      const int &max_infections
				      ){
  // int n_samples = measurement_map_indices.size();
  // int max_infections = cumu_infection_history.size();
  int index;
  double mu = theta["mu"];
  double mu_short = theta["mu_short"];
  //Rcpp::Rcout << "mu_short: " << mu_short << std::endl;
  double senior, inf_map_index, wane;

  for(int i = 0; i < max_infections; ++i){
    //n_inf = cumu_infection_history[i] - 1.0;
    inf_map_index = infection_map_indices[i];
    wane = waning[i];
    senior = seniority[i];
    //if(masked_infection_history[i] > 0){
    for(int k = 0; k < n_samples; ++k){
      index = measurement_map_indices[k] * number_strains + inf_map_index;
      predicted_titres[k] += senior *
	((mu * antigenic_map_long[index]) + 
	 (mu_short * antigenic_map_short[index]) * 
	 wane);    
      //}
    }
  }
}

//[[Rcpp::export]]
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
  
  
  double n_inf, inf_map_index;

  int max_infections = infection_times.size();
  int n_samples = measurement_map_indices.size();

  for(int i = 0; i < max_infections; ++i){
    circulation_time = infection_times[i];
    n_inf = cumu_infection_history[i] - 1.0;
    inf_map_index = infection_map_indices[i];

    if(masked_infection_history[i] > 0){
      for(int ii = i - 1; ii >= 0; --ii){
	if(masked_infection_history[ii] > 0){
	  long_boost = MAX(0, 1.0 - tau*(cumu_infection_history[ii] - 1.0)) * // Antigenic seniority
	    (mu * antigenic_map_long[inf_map_index * 
				     number_strains + infection_map_indices[ii]]);      
      
	  // Short term cross reactive boost
	  short_boost =  MAX(0, 1.0 - tau*(cumu_infection_history[ii] - 1.0)) * // Antigenic seniority
	    (mu_short * antigenic_map_short[inf_map_index * 
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
	long_boost = MAX(0, 1.0 - tau * n_inf) *
	  (mu * antigenic_map_long[measurement_map_indices[k] * 
				   number_strains + inf_map_index]);
    
	// Short term cross reactive boost
	short_boost = MAX(0, 1.0 - tau * n_inf) *
	  (mu_short * antigenic_map_short[measurement_map_indices[k] * 
					  number_strains + inf_map_index]);
    
    
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
				   const NumericVector &seniority,
				   const int &number_strains,
				   const int &n_samples,
				   const int &max_infections,
				   const bool &titre_dependent_boosting,
				   const int &DOB,
				   const Nullable<List> &additional_arguments
				   ){ 
  if (titre_dependent_boosting) {
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
  } else if (additional_arguments.isNotNull()) {
    List _additional_arguments(additional_arguments);
    multiple_infection_strain_dependent(predicted_titres,
					theta, 
					cumu_infection_history,
					//masked_infection_history,
					infection_map_indices,
					measurement_map_indices,
					antigenic_map_long,
					antigenic_map_short,
					waning,
					number_strains,
					_additional_arguments);	
  } else {
    multiple_infection_base_boosting(predicted_titres,
				     theta, 
				     cumu_infection_history,
				     //masked_infection_history,
				     infection_map_indices,
				     measurement_map_indices,
				     antigenic_map_long,
				     antigenic_map_short,
				     waning,
				     seniority,
				     number_strains,
				     n_samples,
				     max_infections);							       
  }
}

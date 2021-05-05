#include <cmath>
#include "wane_function.h"
#include "boosting_functions_fast.h"
#include "helpers.h"

//' Overall model function, fast implementation
//'
//' @param theta NumericVector, the named vector of model parameters
//' @param infection_history_mat IntegerMatrix, the matrix of 1s and 0s showing presence/absence of infection for each possible time for each individual. 
//' @param circulation_times NumericVector, the actual times of circulation that the infection history vector corresponds to
//' @param circulation_times_indices IntegerVector, which entry in the melted antigenic map that these infection times correspond to
//' @param sample_times NumericVector, the times that each blood sample was taken
//' @param rows_per_indiv_in_samples IntegerVector, one entry for each individual. Each entry dictates how many indices through sample_times to iterate per individual (ie. how many sample times does each individual have?)
//' @param cum_nrows_per_individual_in_data IntegerVector, How many cumulative rows in the titre data correspond to each individual? 
//' @param nrows_per_blood_sample IntegerVector, one entry per sample taken. Dictates how many entries to iterate through cum_nrows_per_individual_in_data for each sampling time considered
//' @param measurement_strain_indices IntegerVector, the indices of all measured strains in the melted antigenic map, with one entry per measured titre
//' @param antigenic_map_long NumericVector, the collapsed cross reactivity map for long term boosting, after multiplying by sigma1 see \code{\link{create_cross_reactivity_vector}}
//' @param antigenic_map_short NumericVector, the collapsed cross reactivity map for short term boosting, after multiplying by sigma2, see \code{\link{create_cross_reactivity_vector}}
//' @param antigenic_distances NumericVector, the collapsed cross reactivity map giving euclidean antigenic distances, see \code{\link{create_cross_reactivity_vector}}
//' @param mus NumericVector, if length is greater than one, assumes that strain-specific boosting is used rather than a single boosting parameter
//' @param boosting_vec_indices IntegerVector, same length as circulation_times, giving the index in the vector \code{mus} that each entry should use as its boosting parameter.
//' @param boost_before_infection bool to indicate if calculated titre for that time should be before the infection has occurred, used to calculate titre-mediated immunity
//' @return NumericVector of predicted titres for each entry in measurement_strain_indices
//' @export
//' @family titre_model
// [[Rcpp::export(rng = false)]]
NumericVector titre_data_fast(const NumericVector &theta, 
			      const IntegerMatrix &infection_history_mat, 
			      const NumericVector &circulation_times,
			      const IntegerVector &circulation_times_indices,
			      const NumericVector &sample_times,
			      const IntegerVector &rows_per_indiv_in_samples, // How many rows in titre data correspond to each individual, sample and repeat?
			      const IntegerVector &cum_nrows_per_individual_in_data, // How many rows in the titre data correspond to each individual?
			      const IntegerVector &nrows_per_blood_sample, // Split the sample times and runs for each individual
			      const IntegerVector &measurement_strain_indices, // For each titre measurement, corresponding entry in antigenic map
			      const NumericVector &antigenic_map_long,
			      const NumericVector &antigenic_map_short,
			      const NumericVector &antigenic_distances,	// Currently not doing anything, but has uses for model extensions		      
			      const NumericVector &mus,
			      const IntegerVector &boosting_vec_indices,
			      bool boost_before_infection = false
			      ){
  // Dimensions of structures
  int n = infection_history_mat.nrow();
  int number_strains = infection_history_mat.ncol();
  int total_titres = measurement_strain_indices.size();
  
  // To track how far through the larger vectors we move for each individual
  int index_in_samples;
  int end_index_in_samples;
  int start_index_in_data;
  
  // Only use the infections that actually happened
  IntegerVector infection_history(number_strains);
  LogicalVector indices;
  
  NumericVector infection_times;
  IntegerVector infection_strain_indices_tmp;

  // ====================================================== //
  // =============== SETUP MODEL PARAMETERS =============== //
  // ====================================================== //
  // 1. Extract general parameters that apply to all models
  // Pull out model parameters so only need to allocate once
  double mu = theta["mu"];
  double mu_short = theta["mu_short"];
  double wane = theta["wane"];
  double tau = theta["tau"];
  double min_titre = 0; //theta["min_titre"];

  // 2. Extract model parameters that are for specific mechanisms
  //    set a boolean flag to choose between model versions
  
  // Alternative waning function
  int wane_type = theta["wane_type"]; 
  bool alternative_wane_func = wane_type == 1;
  double kappa;
  double t_change;
 if (alternative_wane_func){
    kappa = theta["kappa"];
    t_change = theta["t_change"];
  }
 
  // Titre dependent boosting
  bool titre_dependent_boosting = theta["titre_dependent"] == 1;
  double gradient;
  double boost_limit;
  if (titre_dependent_boosting) {
    gradient = theta["gradient"];
    boost_limit = theta["boost_limit"];
  }

  // Strain-specific boosting
  bool strain_dep_boost = false;
  if (mus.size() > 1) {
    strain_dep_boost = true;    
  }
  
  // 3. If not using one of the specific mechanism functions, set the base_function flag to TRUE
  bool base_function = !(alternative_wane_func ||
			 titre_dependent_boosting ||
			 strain_dep_boost);

  // To store calculated titres
  NumericVector predicted_titres(total_titres, min_titre);
  // For each individual
  for (int i = 1; i <= n; ++i) {
    infection_history = infection_history_mat(i-1,_);
    indices = infection_history > 0;
    infection_times = circulation_times[indices];
    // Only solve is this individual has had infections
    if (infection_times.size() > 0) {
      infection_strain_indices_tmp = circulation_times_indices[indices];
    
      index_in_samples = rows_per_indiv_in_samples[i-1];
      end_index_in_samples = rows_per_indiv_in_samples[i] - 1;
      start_index_in_data = cum_nrows_per_individual_in_data[i-1];

      // ====================================================== //
      // =============== CHOOSE MODEL TO SOLVE =============== //
      // ====================================================== //
      // Go to sub function - this is where we have options for different models
      // Note, these are in "boosting_functions.cpp"
      if (base_function) {
	titre_data_fast_individual_base(predicted_titres, mu, mu_short,
					wane, tau,
					infection_times,
					infection_strain_indices_tmp,
					measurement_strain_indices,
					sample_times,
					index_in_samples,
					end_index_in_samples,
					start_index_in_data,
					nrows_per_blood_sample,
					number_strains,
					antigenic_map_short,
					antigenic_map_long,
					boost_before_infection);
      } else if (titre_dependent_boosting) {
	titre_data_fast_individual_titredep(predicted_titres, mu, mu_short,
					    wane, tau,
					    gradient, boost_limit,
					    infection_times,
					    infection_strain_indices_tmp,
					    measurement_strain_indices,
					    sample_times,
					    index_in_samples,
					    end_index_in_samples,
					    start_index_in_data,
					    nrows_per_blood_sample,
					    number_strains,
					    antigenic_map_short,
					    antigenic_map_long,
					    boost_before_infection);	
      } else if (strain_dep_boost) {
	titre_data_fast_individual_strain_dependent(predicted_titres, 
						    mus, boosting_vec_indices, 
						    mu_short,
						    wane, tau,
						    infection_times,
						    infection_strain_indices_tmp,
						    measurement_strain_indices,
						    sample_times,
						    index_in_samples,
						    end_index_in_samples,
						    start_index_in_data,
						    nrows_per_blood_sample,
						    number_strains,
						    antigenic_map_short,
						    antigenic_map_long,
						    boost_before_infection);
      } else if(alternative_wane_func) {
	titre_data_fast_individual_wane2(predicted_titres, mu, mu_short,
					 wane, tau,
					 kappa, t_change,
					 infection_times,
					 infection_strain_indices_tmp,
					 measurement_strain_indices,
					 sample_times,
					 index_in_samples,
					 end_index_in_samples,
					 start_index_in_data,
					 nrows_per_blood_sample,
					 number_strains,
					 antigenic_map_short,
					 antigenic_map_long,
					 boost_before_infection);
      } else {
	titre_data_fast_individual_base(predicted_titres, mu, mu_short,
					wane, tau,
					infection_times,
					infection_strain_indices_tmp,
					measurement_strain_indices,
					sample_times,
					index_in_samples,
					end_index_in_samples,
					start_index_in_data,
					nrows_per_blood_sample,
					number_strains,
					antigenic_map_short,
					antigenic_map_long,
					boost_before_infection);
      }
     
    }
  }
  return(predicted_titres);
}

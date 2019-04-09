#include <cmath>
#include "wane_function.h"
#include "boosting_functions_fast.h"
#include "helpers.h"

//' Overall model function, fast implementation
//'
//' See documentation for \code{\link{titre_data_group}}, as the interface is almost identical
//' @inheritParams titre_data_group
//' @param mus NumericVector, if length is greater than one, assumes that strain-specific boosting is used rather than a single boosting parameter
//' @param boosting_vec_indices IntegerVector, same length as circulation_times, giving the index in the vector \code{mus} that each entry should use as its boosting parameter.
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
			      const NumericVector &mus,
			      const IntegerVector &boosting_vec_indices
			      ){
  // Dimensions of structures
  int n = infection_history_mat.nrow();
  int number_strains = infection_history_mat.ncol();
  int total_titres = measurement_strain_indices.size();
  int max_infections;
  int n_titres;
  
  // To track how far through the larger vectors we move for each individual
  int index_in_samples;
  int end_index_in_samples;
  int number_samples;
  int start_index_in_data;
  int end_index_in_data;
  int tmp_titre_index;
  int inf_map_index;
  int index;

  // 
  double sampling_time;
  double time;
  double n_inf;
  
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
  NumericVector predicted_titres(total_titres);
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
      number_samples = end_index_in_samples - index_in_samples;      
      start_index_in_data = cum_nrows_per_individual_in_data[i-1];

      // ====================================================== //
      // =============== CHOOSE MODEL TO SOLVE =============== //
      // ====================================================== //
      // Go to sub function - this is where we'd have options for different models
      // Note, this is in "boosting_functions.cpp"
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
					antigenic_map_long);
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
					    antigenic_map_long);	
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
						    antigenic_map_long);
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
					 antigenic_map_long);
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
					antigenic_map_long);
      }
     
    }
  }
  return(predicted_titres);
}

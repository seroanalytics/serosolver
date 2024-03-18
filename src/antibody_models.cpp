#include <cmath>
#include "antibody_models_individual.h"
#include "helpers.h"

//' Overall model function, fast implementation
//'
//' @param theta NumericVector, the named vector of model parameters
//' @param infection_history_mat IntegerMatrix, the matrix of 1s and 0s showing presence/absence of infection for each possible time for each individual. 
//' @param possible_exposure_times NumericVector, the time periods that the infection history vector corresponds to
//' @param possible_exposure_times_indices IntegerVector, which entry in the melted antigenic map that each possible infection time corresponds to
//' @param sample_times NumericVector, the times that each blood sample was taken
//' @param type_data_start IntegerVector, one entry for each unique individual. Each entry gives the starting index for each individual of the data frame `unique(antibody_data[,c("individual","biomarker_group")])`.
//' @param biomarker_groups IntegerVector, result of `unique(antibody_data[,c("individual","biomarker_group")])$biomarker_group`
//' @param sample_data_start IntegerVector, one entry for each unique individual and observation type combination. Each entry dictates how many indices through sample_times to iterate per individual and observation type (ie. how many sample times does each individual have?)
//' @param antibody_data_start IntegerVector, How many cumulative rows in the antibody data correspond to each unique individual and observation type combination? 
//' @param nrows_per_sample IntegerVector, one entry per sample taken. Dictates how many entries to iterate through cum_nrows_per_individual_in_data for each sampling time considered
//' @param biomarker_id_indices IntegerVector, the indices of all measured biomarkers in the melted antigenic map, with one entry per measured biomarker
//' @param antigenic_map_long arma::mat, the collapsed cross reactivity map for long term boosting, after multiplying by cr_long see \code{\link{create_cross_reactivity_vector}}
//' @param antigenic_map_short arma::mat, the collapsed cross reactivity map for short term boosting, after multiplying by cr_short, see \code{\link{create_cross_reactivity_vector}}
//' @param antigenic_distances NumericVector, the collapsed cross reactivity map giving euclidean antigenic distances, see \code{\link{create_cross_reactivity_vector}}
//' @param boost_before_infection bool to indicate if calculated antibody level for that time should be before the infection has occurred, used to calculate antibody-mediated immunity
//' @return NumericVector of predicted antibody levels for each entry in biomarker_id_indices
//' @export
//' @family antibody_models
// [[Rcpp::export(rng = false)]]
NumericVector antibody_model(const NumericVector &theta, 
                              const IntegerVector &unique_theta_indices,
                              const IntegerVector &unique_biomarker_groups,
			      const IntegerMatrix &infection_history_mat, 
			      const NumericVector &possible_exposure_times,
			      const IntegerVector &possible_exposure_times_indices,
			      const NumericVector &sample_times,
			      
			      const IntegerVector &type_data_start,
			      const IntegerVector &biomarker_groups,
			      const IntegerVector &sample_data_start, 
			      const IntegerVector &antibody_data_start, 
			      const IntegerVector &nrows_per_sample, // Split the sample times for each individual
			      
			      const IntegerVector &biomarker_id_indices, // For each antibody measurement, corresponding entry in antigenic map
			      const arma::mat &antigenic_map_long,
			      const arma::mat &antigenic_map_short,
			      const NumericVector &antigenic_distances,	// Currently not doing anything, but has uses for model extensions
			      const Nullable<NumericVector> &starting_antibody_levels = R_NilValue,
			      bool boost_before_infection = false
			      ){
  // Dimensions of structures
  int n = infection_history_mat.nrow();
  int number_possible_infections = infection_history_mat.ncol();
  int n_measurements = biomarker_id_indices.size();

  // To track how far through the larger vectors we move for each individual
  int biomarker_group=0;
  int type_start;
  int type_end;
  int start_index_in_samples;
  int end_index_in_samples;
  int start_index_in_data;
  
  // Only use the infections that actually happened
  //IntegerVector infection_history(number_possible_infections);
  LogicalVector indices;
  
  NumericVector infection_times;
  IntegerVector infection_times_indices_tmp;

  // ====================================================== //
  // =============== SETUP MODEL PARAMETERS =============== //
  // ====================================================== //
  // 1. Extract general parameters that apply to all models
  // Pull out model parameters so only need to allocate once
  int n_types = unique_biomarker_groups.size();
  int n_theta = unique_theta_indices.size();
  
  // Base model parameters
  NumericVector boost_long_parameters(n_types);
  NumericVector boost_short_parameters(n_types);
  NumericVector wane_short_short_parameters(n_types);
  NumericVector obs_sd_parameters(n_types);
  NumericVector antigenic_seniority_parameters(n_types);
  
  //NumericVector max_measurements(n_types);
  //NumericVector min_measurements(n_types);
  double min_measurement = 0;
  
  int boost_long_index = unique_theta_indices("boost_long");
  int boost_short_index = unique_theta_indices("boost_short");
  int wane_short_index = unique_theta_indices("wane_short");
  int antigenic_seniority_index = unique_theta_indices("antigenic_seniority");
  int error_index = unique_theta_indices("obs_sd");
  
  //NumericVector antigenic_map_long_tmp(antigenic_map_long.nrow());
  //NumericVector antigenic_map_short_tmp(antigenic_map_short.nrow());
  
  //int min_index = unique_theta_indices("min_measurement");
  //int max_index = unique_theta_indices("max_measurement");
  
  // Titre-dependent boosting function
  IntegerVector antibody_dependent_boosting(n_types);
  NumericVector gradients(n_types); 
  NumericVector boost_limits(n_types);   
  
  int antibody_dependent_boosting_index = unique_theta_indices("antibody_dependent_boosting");
  int gradient_index = -1;
  int boost_limit_index = -1;

  // Create vectors of model parameters for each of the observation types

  for(int x = 0; x < n_types; ++x){

      boost_long_parameters(x) = theta(boost_long_index + x*n_theta);
      boost_short_parameters(x) = theta(boost_short_index + x*n_theta);
      wane_short_short_parameters(x) = theta(wane_short_index + x*n_theta);
      antigenic_seniority_parameters(x) = theta(antigenic_seniority_index + x*n_theta);
      obs_sd_parameters(x) = theta(error_index + x*n_theta);
      
      //min_measurements(x) = theta(min_index + x*n_theta);
      //max_measurements(x) = theta(max_index + x*n_theta);
      
      // Titre dependent boosting
      antibody_dependent_boosting(x) = theta(antibody_dependent_boosting_index+ x*n_theta);
      
      if(antibody_dependent_boosting(x) == 1) {
          gradient_index = unique_theta_indices("gradient");
          boost_limit_index = unique_theta_indices("boost_limit");  
          gradients(x) = theta(gradient_index + x*n_theta);
          boost_limits(x) = theta(boost_limit_index + x*n_theta);
      }
  }

  NumericVector predicted_antibody_levels = get_starting_antibody_levels(n_measurements,min_measurement,starting_antibody_levels);
 // antigenic_map_long_tmp = antigenic_map_long(_,biomarker_group);
  //antigenic_map_short_tmp = antigenic_map_short(_,biomarker_group);
  
  // For each individual
  for (int i = 1; i <= n; ++i) {
    indices = infection_history_mat(i-1,_) > 0;
    
    infection_times = possible_exposure_times[indices];
    // Only solve is this individual has had infections
    if (infection_times.size() > 0) {
      infection_times_indices_tmp = possible_exposure_times_indices[indices];
    
      // Start end end location of the type_data matrix
      type_start = type_data_start(i-1);
      type_end = type_data_start(i)-1;
    
      // For each observation type solved for this individual
          for(int index = type_start; index <= type_end; ++index){
              biomarker_group = biomarker_groups(index)-1;

            start_index_in_samples = sample_data_start(index);
            end_index_in_samples = sample_data_start(index+1) - 1;
            start_index_in_data = antibody_data_start(start_index_in_samples);
    
      
    
    
            // ====================================================== //
            // =============== CHOOSE MODEL TO SOLVE =============== //
            // ====================================================== //
            // Go to sub function - this is where we have options for different models
            // Note, these are in "boosting_functions.cpp"
            if (antibody_dependent_boosting(biomarker_group)) {
              antibody_dependent_boosting_model_individual(
                predicted_antibody_levels, 
                boost_long_parameters(biomarker_group), 
                boost_short_parameters(biomarker_group),
	              wane_short_short_parameters(biomarker_group), 
	              antigenic_seniority_parameters(biomarker_group),
  					    gradients(biomarker_group), 
  					    boost_limits(biomarker_group),
  					    infection_times,
  					    infection_times_indices_tmp,
  					    biomarker_id_indices,
  					    sample_times,
  					    start_index_in_samples,
  					    end_index_in_samples,
  					    start_index_in_data,
  					    nrows_per_sample,
  					    number_possible_infections,
  					    antigenic_map_short.colptr(biomarker_group),
  					    antigenic_map_long.colptr(biomarker_group),
  					    boost_before_infection);	
            } else {
              antibody_data_model_individual(
            	        predicted_antibody_levels, 
                      boost_long_parameters(biomarker_group), 
                      boost_short_parameters(biomarker_group),
            					wane_short_short_parameters(biomarker_group), 
            					antigenic_seniority_parameters(biomarker_group),
            					infection_times,
            					infection_times_indices_tmp,
            					biomarker_id_indices,
            					sample_times,
            					start_index_in_samples,
            					end_index_in_samples,
            					start_index_in_data,
            					nrows_per_sample,
            					number_possible_infections,
            					antigenic_map_short.colptr(biomarker_group),
            					antigenic_map_long.colptr(biomarker_group),
            					boost_before_infection);
              }
          }
      }
  }
  return(predicted_antibody_levels);
}

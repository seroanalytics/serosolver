#include <cmath>
#include "antibody_models_individual.h"
#include "helpers.h"

//' Overall model function, fast implementation
//'
//' Overall model function, fast implementation
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
NumericVector antibody_model(const NumericMatrix theta, 
                             
                              const IntegerVector &unique_theta_indices,
                              const IntegerVector &unique_biomarker_groups,
			      const IntegerMatrix &infection_history_mat, 
			      const IntegerVector &infection_history_mat_indices,
			      
			      const IntegerVector &indiv_theta_groups,
			      
			      const NumericVector &possible_exposure_times,
			      const IntegerVector &possible_exposure_times_indices,
			      const NumericVector &sample_times,
			      
			      const IntegerVector &type_data_start,
			      const IntegerVector &biomarker_groups,
			      const IntegerVector &sample_data_start, 
			      const IntegerVector &antibody_data_start, 
			      const IntegerVector &nrows_per_sample, // Split the sample times for each individual
			      
			      const IntegerVector &biomarker_id_indices, // For each antibody measurement, corresponding entry in antigenic map
			      const IntegerVector &start_level_indices, // For each individual/biomarker group/biomarker id combo, we need to know the starting antibody level
			      const NumericVector &starting_antibody_levels,
			      const NumericVector &births,
			      
			      const arma::cube &antigenic_map_long,
			      const arma::cube &antigenic_map_short,
			      const NumericVector &antigenic_distances,	// Currently not doing anything, but has uses for model extensions
			      const bool timevarying_groups = false,
			      bool boost_before_infection = false
			      ){
  // Dimensions of structures
  int n = infection_history_mat.nrow();
  int number_possible_exposures = possible_exposure_times.size();
  int n_measurements = biomarker_id_indices.size();
  //NumericVector predicted_antibody_levels(n_measurements, 0.0);
  //return predicted_antibody_levels;
  
  
  // To track how far through the larger vectors we move for each individual
  int biomarker_group=0;
  
  // Tracking an individual's demographic group
  IntegerVector groups;
  int group = 0;
  int birth_group=0;
  
  int type_start;
  int type_end;
  int start_index_in_samples;
  int end_index_in_samples;
  int start_index_in_data;
  
  // Only use the infections that actually happened
  //IntegerVector infection_history(number_possible_exposures);
  LogicalVector indices;
  
  NumericVector infection_times;
  IntegerVector infection_times_indices_tmp;
  IntegerVector use_indices;

  // ====================================================== //
  // =============== SETUP MODEL PARAMETERS =============== //
  // ====================================================== //
  // 1. Extract general parameters that apply to all models
  // Pull out model parameters so only need to allocate once
  int n_types = unique_biomarker_groups.size();
  int n_theta = unique_theta_indices.size();
  int n_groups = theta.nrow();
  
  // Base model parameters
  // Matrix with entry for each group and biomarker type
  NumericMatrix boost_long_parameters(n_groups, n_types);
  NumericMatrix boost_short_parameters(n_groups, n_types);
  NumericMatrix boost_delay_parameters(n_groups, n_types);
  NumericMatrix wane_short_parameters(n_groups, n_types);
  NumericMatrix wane_long_parameters(n_groups, n_types);
  NumericMatrix wane_maternal_parameters(n_groups, n_types);
  NumericMatrix obs_sd_parameters(n_groups, n_types);
  NumericMatrix antigenic_seniority_parameters(n_groups, n_types);
  NumericMatrix min_measurements(n_groups, n_types);
  
  int boost_long_index = unique_theta_indices("boost_long");
  int boost_short_index = unique_theta_indices("boost_short");
  int boost_delay_index = unique_theta_indices("boost_delay");
  int wane_short_index = unique_theta_indices("wane_short");
  int wane_long_index = unique_theta_indices("wane_long");
  int wane_maternal_index = unique_theta_indices("wane_maternal");
  int antigenic_seniority_index = unique_theta_indices("antigenic_seniority");
  int error_index = unique_theta_indices("obs_sd");
  int min_index = unique_theta_indices("min_measurement");
  
  // Titre-dependent boosting function
  IntegerMatrix antibody_dependent_boosting(n_groups,n_types);
  NumericMatrix gradients(n_groups,n_types); 
  NumericMatrix boost_limits(n_groups,n_types);   
  
  int antibody_dependent_boosting_index = unique_theta_indices("antibody_dependent_boosting");
  int gradient_index = -1;
  int boost_limit_index = -1;

  // Create vectors of model parameters for each of the observation types
  NumericVector predicted_antibody_levels(n_measurements, 0.0);

  for(int g = 0; g < n_groups; ++g){
    for(int x = 0; x < n_types; ++x){
  
        boost_long_parameters(g,x) = theta(g,boost_long_index + x*n_theta);
        boost_short_parameters(g,x) = theta(g,boost_short_index + x*n_theta);
        boost_delay_parameters(g,x) = theta(g,boost_delay_index + x*n_theta);
        wane_short_parameters(g,x) = theta(g,wane_short_index + x*n_theta);
        wane_long_parameters(g,x) = theta(g,wane_long_index + x*n_theta);
        wane_maternal_parameters(g,x) = theta(g,wane_maternal_index + x*n_theta);  // Maternal antibody waning
        antigenic_seniority_parameters(g,x) = theta(g,antigenic_seniority_index + x*n_theta);
        obs_sd_parameters(g,x) = theta(g,error_index + x*n_theta);
        
        min_measurements(g,x) = theta(g,min_index + x*n_theta);
        
        // Titre dependent boosting
        antibody_dependent_boosting(g,x) = theta(g,antibody_dependent_boosting_index+ x*n_theta);
        
        
        if(antibody_dependent_boosting(g,x) == 1) {
            gradient_index = unique_theta_indices("gradient");
            boost_limit_index = unique_theta_indices("boost_limit");  
            gradients(g,x) = theta(g,gradient_index + x*n_theta);
            boost_limits(g,x) = theta(g,boost_limit_index + x*n_theta);
        }
      }
  }

  // For each individual
  for (int i = 1; i <= n; ++i) {
    indices = infection_history_mat(i-1,_) > 0;
    
    use_indices =infection_history_mat_indices[indices];
    infection_times = possible_exposure_times[use_indices];
    infection_times_indices_tmp = possible_exposure_times_indices[use_indices];	 
    
    // If trying to solve how demographic groups change over time, then need to extract the demographic group indices for this individuals
    if(timevarying_groups){
      // Treat the indiv_theta_groups vector as a flattened matrix -- this individual's demography vector is from [i-1,2] to [i-1, number_possible_exposures + 1] (the +1 is for birth group) 
      groups = indiv_theta_groups[Range((number_possible_exposures+1)*(i-1) + 1, (number_possible_exposures+1)*i - 1)];
      //Rcpp::Rcout << "Indiv: " << i << std::endl;
      //Rcpp::Rcout << "Groups: " << groups << std::endl;
      
      // Then subset to groups with infections
      groups = groups[use_indices];
      
      //Rcpp::Rcout << "Groups subset: " << groups << std::endl;
      // Birth group is [i-1, 1]
      birth_group = indiv_theta_groups[(number_possible_exposures+1)*(i-1)];
      //Rcpp::Rcout << "Infection times: " << infection_times << std::endl;
      
      //Rcpp::Rcout << "Birth group: " << birth_group << std::endl;
      //Rcpp::Rcout << "Theta: " << theta << std::endl;
    } else {
      // Otherwise, the individual's group is unchanged
      group = indiv_theta_groups[i-1];
    }
   
    // Only solve is this individual has had infections or if there is long-term waning
    // if (wane_long_parameters(biomarker_group) > 0 || infection_times.size() > 0) {
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
            if (antibody_dependent_boosting(group,biomarker_group)) {
              antibody_dependent_boosting_model_individual(
                predicted_antibody_levels, 
                boost_long_parameters(group,biomarker_group), 
                boost_short_parameters(group,biomarker_group),
	              wane_short_parameters(group,biomarker_group), 
	              antigenic_seniority_parameters(group,biomarker_group),
  					    gradients(group,biomarker_group), 
  					    boost_limits(group,biomarker_group),
  					    infection_times,
  					    infection_times_indices_tmp,
  					    biomarker_id_indices,
  					    sample_times,
  					    start_index_in_samples,
  					    end_index_in_samples,
  					    start_index_in_data,
  					    nrows_per_sample,
  					    number_possible_exposures,
  					    antigenic_map_short.slice(group).colptr(biomarker_group),
  					    antigenic_map_long.slice(group).colptr(biomarker_group),
  					    boost_before_infection);
            } else if(timevarying_groups){
              antibody_data_model_individual_timevarying(
                predicted_antibody_levels, 
                starting_antibody_levels,
                births,
                boost_long_parameters, 
                boost_short_parameters,
                boost_delay_parameters,
                wane_short_parameters, 
                wane_long_parameters, 
                wane_maternal_parameters,
                antigenic_seniority_parameters,
                infection_times,
                groups,
                birth_group,
                infection_times_indices_tmp,
                biomarker_id_indices,
                start_level_indices,
                sample_times,
                start_index_in_samples,
                end_index_in_samples,
                start_index_in_data,
                nrows_per_sample,
                number_possible_exposures,
                antigenic_map_short,
                antigenic_map_long,
                biomarker_group,
                min_measurements,
                boost_before_infection);
              
            } else {
              antibody_data_model_individual_new(
            	        predicted_antibody_levels, 
            	        starting_antibody_levels,
            	        births,
                      boost_long_parameters(group,biomarker_group), 
                      boost_short_parameters(group,biomarker_group),
                      boost_delay_parameters(group,biomarker_group),
            					wane_short_parameters(group,biomarker_group), 
            					wane_long_parameters(group,biomarker_group), 
            					wane_maternal_parameters(group,biomarker_group), 
                      antigenic_seniority_parameters(group,biomarker_group),
            					infection_times,
            					infection_times_indices_tmp,
            					biomarker_id_indices,
            					start_level_indices,
            					sample_times,
            					start_index_in_samples,
            					end_index_in_samples,
            					start_index_in_data,
            					nrows_per_sample,
            					number_possible_exposures,
            					antigenic_map_short.slice(group).colptr(biomarker_group),
            					antigenic_map_long.slice(group).colptr(biomarker_group),
            					boost_before_infection,
            					min_measurements(group,biomarker_group));
              }
         }
  }
  return(predicted_antibody_levels);
  
}
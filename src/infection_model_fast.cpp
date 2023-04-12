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
//' @param type_data_start IntegerVector, one entry for each unique individual. Each entry gives the starting index for each individual of the data frame `unique(titre_dat[,c("individual","obs_type")])`.
//' @param obs_types IntegerVector, result of `unique(titre_dat[,c("individual","obs_type")])$obs_type`
//' @param sample_data_start IntegerVector, one entry for each unique individual and observation type combination. Each entry dictates how many indices through sample_times to iterate per individual and observation type (ie. how many sample times does each individual have?)
//' @param titre_data_start IntegerVector, How many cumulative rows in the titre data correspond to each unique individual and observation type combination? 
//' @param nrows_per_sample IntegerVector, one entry per sample taken. Dictates how many entries to iterate through cum_nrows_per_individual_in_data for each sampling time considered
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
                              const IntegerVector &unique_theta_indices,
                              const IntegerVector &unique_obs_types,
			      const IntegerMatrix &infection_history_mat, 
			      const NumericVector &circulation_times,
			      const IntegerVector &circulation_times_indices,
			      const NumericVector &sample_times,
			      
			      const IntegerVector &type_data_start,
			      const IntegerVector &obs_types,
			      const IntegerVector &sample_data_start, 
			      const IntegerVector &titre_data_start, 
			      const IntegerVector &nrows_per_sample, // Split the sample times for each individual
			      
			      const IntegerVector &measurement_strain_indices, // For each titre measurement, corresponding entry in antigenic map
			      const NumericMatrix &antigenic_map_long,
			      const NumericMatrix &antigenic_map_short,
			      const NumericVector &antigenic_distances,	// Currently not doing anything, but has uses for model extensions		      
			      const NumericVector &mus_strain_dep,
			      const IntegerVector &boosting_vec_indices,
			      bool boost_before_infection = false
			      ){
  // Dimensions of structures
  int n = infection_history_mat.nrow();
  int number_strains = infection_history_mat.ncol();
  int total_titres = measurement_strain_indices.size();
  
  // To track how far through the larger vectors we move for each individual
  int obs_type=0;
  int type_start;
  int type_end;
  int start_index_in_samples;
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
  int n_types = unique_obs_types.size();
  int n_theta = unique_theta_indices.size();
  
  // Base model parameters
  NumericVector mus(n_types);
  NumericVector mu_shorts(n_types);
  NumericVector wanes(n_types);
  NumericVector errors(n_types);
  NumericVector taus(n_types);
  
  NumericVector min_titres(n_types);
  double min_titre = 0;
  
  int mu_index = unique_theta_indices["mu"];
  int mu_short_index = unique_theta_indices["mu_short"];
  int wane_index = unique_theta_indices["wane"];
  int tau_index = unique_theta_indices["tau"];
  int error_index = unique_theta_indices["error"];

  // Titre-dependent boosting function
  bool titre_dependent_boosting=false;
  NumericVector gradients(n_types); 
  NumericVector boost_limits(n_types);   
  
  int titre_dependent_boosting_index = unique_theta_indices["titre_dependent"];
  int gradient_index = unique_theta_indices["gradient"];
  int boost_limit_index = unique_theta_indices["boost_limit"];  
  
  // Strain-specific boosting
  bool strain_dep_boost = false;
  if (mus_strain_dep.size() > 1) {
      strain_dep_boost = true;    
  }
  
  // Create vectors of model parameters for each of the observation types
  for(int x = 0; x < n_types; ++x){
      mus[x] = theta[mu_index + x*n_theta];
      mu_shorts[x] = theta[mu_short_index + x*n_theta];
      wanes[x] = theta[wane_index + x*n_theta];
      taus[x] = theta[tau_index + x*n_theta];
      errors[x] = theta[error_index + x*n_theta];
      
      // Titre dependent boosting
      //titre_dependent_boosting = theta[titre_dependent_boosting_index+ x*n_theta];
     // if (titre_dependent_boosting) {
       //   gradients(x) = theta[gradient_index + x*n_theta];
       //   boost_limits(x) = theta[boost_limit_index + x*n_theta];
      //}
  }

  // To store calculated titres
  NumericVector predicted_titres(total_titres, min_titre);
  //Rcpp::Rcout << "Sample times: " << sample_times << std::endl;
  // For each individual
  for (int i = 1; i <= n; ++i) {
      //Rcpp::Rcout << std::endl << "Individual: " << i << std::endl;
    infection_history = infection_history_mat(i-1,_);
    indices = infection_history > 0;
    //Rcpp::Rcout << "Indices: " << indices << std::endl;
    
    infection_times = circulation_times[indices];
    // Only solve is this individual has had infections
    if (infection_times.size() > 0) {
        //Rcpp::Rcout << std::endl << "Infection times: "  << infection_times << std::endl;
        
      infection_strain_indices_tmp = circulation_times_indices[indices];
    
    // Start end end location of the type_data matrix
    type_start = type_data_start[i-1];
    type_end = type_data_start[i]-1;
    
    //Rcpp::Rcout << "Type start: " << type_start << std::endl;
    //Rcpp::Rcout << "Type end: " << type_end << std::endl << std::endl;
    
    // For each observation type solved for this individual
        for(int index = type_start; index <= type_end; ++index){
            obs_type = obs_types[index]-1;
           //Rcpp::Rcout << "Observation type: " << obs_type << std::endl;
            
          start_index_in_samples = sample_data_start[index];
          end_index_in_samples = sample_data_start[index+1] - 1;
          start_index_in_data = titre_data_start[start_index_in_samples];
    /*
    Rcpp::Rcout << "Start index in samples: " << start_index_in_samples << std::endl;
    Rcpp::Rcout << "End index in samples: " << end_index_in_samples << std::endl;
    Rcpp::Rcout << "Start index in data: " << start_index_in_data << std::endl;
    Rcpp::Rcout << "Nrows in samples: " << nrows_per_sample << std::endl;
    
    Rcpp::Rcout << "Mu: " << mus(obs_type) << std::endl;
    Rcpp::Rcout << "Mu short: " << mu_shorts(obs_type) << std::endl;
    */
          // ====================================================== //
          // =============== CHOOSE MODEL TO SOLVE =============== //
          // ====================================================== //
          // Go to sub function - this is where we have options for different models
          // Note, these are in "boosting_functions.cpp"
          if (titre_dependent_boosting) {
              titre_data_fast_individual_titredep(predicted_titres, 
                                mus(obs_type), 
                                mu_shorts(obs_type),
				                wanes(obs_type), 
				                taus(obs_type),
        					    gradients(obs_type), 
        					    boost_limits(obs_type),
        					    infection_times,
        					    infection_strain_indices_tmp,
        					    measurement_strain_indices,
        					    sample_times,
        					    start_index_in_samples,
        					    end_index_in_samples,
        					    start_index_in_data,
        					    nrows_per_sample,
        					    number_strains,
        					    antigenic_map_short(_,obs_type),
        					    antigenic_map_long(_,obs_type),
        					    boost_before_infection);	
          } else if (strain_dep_boost) {
        	titre_data_fast_individual_strain_dependent(predicted_titres, 
        						    mus(obs_type), 
        						    boosting_vec_indices, 
        						    mu_shorts(obs_type),
        						    wanes(obs_type), 
        						    taus(obs_type),
        						    infection_times,
        						    infection_strain_indices_tmp,
        						    measurement_strain_indices,
        						    sample_times,
        						    start_index_in_samples,
        						    end_index_in_samples,
        						    start_index_in_data,
        						    nrows_per_sample,
        						    number_strains,
        						    antigenic_map_short(_,obs_type),
        						    antigenic_map_long(_,obs_type),
        						    boost_before_infection);
          } else {
        	titre_data_fast_individual_base(predicted_titres, 
                            mus(obs_type), 
                            mu_shorts(obs_type),
        					wanes(obs_type), 
        					taus(obs_type),
        					infection_times,
        					infection_strain_indices_tmp,
        					measurement_strain_indices,
        					sample_times,
        					start_index_in_samples,
        					end_index_in_samples,
        					start_index_in_data,
        					nrows_per_sample,
        					number_strains,
        					antigenic_map_short(_,obs_type),
        					antigenic_map_long(_,obs_type),
        					boost_before_infection);
          }
        }
    }
  }
  return(predicted_titres);
}

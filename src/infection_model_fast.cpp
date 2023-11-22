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
  double sigma_s = theta["sigma2"];
  double sigma_l = theta["sigma1"];

  // Pre-compute waning and antigenic seniority contributions for speed
  NumericVector wane_amounts(number_strains);
  NumericVector seniority_amounts(number_strains);
  double t = 0;
  for(int index = 0; index < number_strains; ++index){
      wane_amounts[index] = MAX(0, 1.0 - (wane*t));
      seniority_amounts[index] = MAX(0, 1.0 - (tau*t));
      t = t + 1.0;
  }
  
  // 2. Extract model parameters that are for specific mechanisms
  //    set a boolean flag to choose between model versions
  // Titre dependent boosting
  bool titre_dependent_boosting = theta["titre_dependent"] == 1;
  double gradient;
  double boost_limit;
  if (titre_dependent_boosting) {
    gradient = theta["gradient"];
    boost_limit = theta["boost_limit"];
  }
  
  // Complex cross-reactivity model
  bool complex_cr = theta["complex_cr"] == 1;
    double sigma_birth_mod_s;
    double sigma_birth_mod_l;
    double sigma_future_mod_s;
    double sigma_future_mod_l;
    if(complex_cr){
        sigma_birth_mod_s = theta["sigma_birth_mod_s"];
        sigma_birth_mod_l= theta["sigma_birth_mod_l"];
        sigma_future_mod_s= theta["sigma_future_mod_s"];
        sigma_future_mod_l= theta["sigma_future_mod_l"];
    }
    
  // 3. If not using one of the specific mechanism functions, set the base_function flag to TRUE
  bool base_function = !(complex_cr ||  titre_dependent_boosting);

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
        					wane_amounts, seniority_amounts,
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
      } else if(complex_cr) {
          titre_data_fast_individual_complex_cr(predicted_titres, mu, mu_short,
                                          wane, tau,sigma_s, sigma_l,
                                          sigma_birth_mod_s, sigma_birth_mod_l,
                                          sigma_future_mod_s, sigma_future_mod_l,
                                          wane_amounts, seniority_amounts,
                                          infection_times,
                                          infection_strain_indices_tmp,
                                          measurement_strain_indices,
                                          sample_times,
                                          index_in_samples,
                                          end_index_in_samples,
                                          start_index_in_data,
                                          nrows_per_blood_sample,
                                          number_strains,
                                          antigenic_distances,
                                          boost_before_infection);
      } else {
        	titre_data_fast_individual_base(predicted_titres, mu, mu_short,
        					wane, tau,
        					wane_amounts, seniority_amounts,
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

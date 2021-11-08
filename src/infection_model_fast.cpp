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
NumericVector titre_data_fast(NumericVector theta,
			      const IntegerMatrix infection_history_mat,
            const DataFrame vaccination_histories,
			      const NumericVector circulation_times,
			      const IntegerVector circulation_times_indices,
			      const NumericVector sample_times,
			      const IntegerVector rows_per_indiv_in_samples, // How many rows in titre data correspond to each individual, sample and repeat?
			      const IntegerVector cum_nrows_per_individual_in_data, // How many rows in the titre data correspond to each individual?
			      const IntegerVector nrows_per_blood_sample, // Split the sample times and runs for each individual
			      const IntegerVector measurement_strain_indices, // For each titre measurement, corresponding entry in antigenic map
			      const NumericVector antigenic_map_long, // add in antigenic_map_long_vac ?
			      const NumericVector antigenic_map_short,
            const NumericVector antigenic_map_long_vac, 
			      const NumericVector antigenic_map_short_vac,
			      const NumericVector antigenic_distances,	// Currently not doing anything, but has uses for model extensions		      
			      const NumericVector mus,
			      const IntegerVector boosting_vec_indices,
			      bool boost_before_infection = false
			      ){
  // Dimensions of structures
  //Rcpp::Rcout << "titre_data_fast: START" << std::endl;
  int n = infection_history_mat.nrow();
  int number_strains = infection_history_mat.ncol();
  int total_titres = measurement_strain_indices.size();
   // Rcpp::Rcout << "Here 1" << std::endl;

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
//  Rcpp::Rcout << "Sec 1" << std::endl;

  double mu = theta["mu"];
  double mu_short = theta["mu_short"];
  double wane = theta["wane"];
  double tau = theta["tau"];
  double min_titre = 0; //theta["min_titre"];
    //Rcpp::Rcout << "Here 2" << std::endl;

  // 2. Extract model parameters that are for specific mechanisms
  //    set a boolean flag to choose between model versions
 // Rcpp::Rcout << "Sec 2" << std::endl;

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

  // Vaccination stuff  
  int vac_flag_int = theta["vac_flag"];
  bool vac_flag = vac_flag_int == 1;
  double mu_vac = 0;
  double mu_short_vac = 0;
  double wane_vac = 0;
  double tau_prev_vac = 0;
  if (vac_flag) {
        mu_vac = theta["mu_vac"];
        mu_short_vac = theta["mu_short_vac"];
        wane_vac = theta["wane_vac"];
        tau_prev_vac = theta["tau_prev_vac"];
  }


  // 3. If not using one of the specific mechanism functions, set the base_function flag to TRUE
  //Rcpp::Rcout << "Sec 3" << std::endl;
  //    Rcpp::Rcout << "Sec 3" << std::endl;
  bool base_function = !(alternative_wane_func ||
			 titre_dependent_boosting ||
			 strain_dep_boost);
  
    LogicalVector vac_flag_ind;
    NumericVector vaccination_times;
    IntegerVector vaccination_strain_indices_tmp;
    NumericVector vaccinations_previous;
      
    std::vector<bool> vac_flag_ind_cpp;
    std::vector<int> vaccination_times_cpp;
    std::vector<int> vaccinations_previous_cpp;
    std::vector<int> vaccination_strain_indices_tmp_cpp;

    int Ni = vaccination_histories.nrows();

    NumericVector individuals_vacc_vec(Ni);
    NumericVector vac_virus_vec(Ni);
    NumericVector vac_time_vec(Ni);
    NumericVector vac_flag_vec(Ni);
    NumericVector prev_vac_vec(Ni);

      

    if (vac_flag) { 
      individuals_vacc_vec = vaccination_histories[0];
      vac_virus_vec = vaccination_histories[1];
      vac_time_vec = vaccination_histories[2];
      vac_flag_vec = vaccination_histories[3];
      prev_vac_vec = vaccination_histories[4];
    } else {

    }

  std::vector<bool> indices_vac_individ(Ni);
    // change these to cpp values
  std::vector<double > vac_virus_vec_ind;
  std::vector<double > vac_time_vec_ind;
  std::vector<double > vaccination_strains;
  std::vector<double > prev_vac_ind;    

  // To store calculated titres
  LogicalVector vac_flag_ind_bool;
  NumericVector predicted_titres(total_titres, min_titre);
  // For each individual (n number of individuals)
  bool indiv_indic;
  for (int i = 1; i <= n; ++i) {
   // Rcpp::Rcout << "Number: " << i << std::endl;

    vac_virus_vec_ind.clear();
    vac_time_vec_ind.clear();
    vac_flag_ind_cpp.clear();
    prev_vac_ind.clear();
    vaccination_times_cpp.clear();
    vaccinations_previous_cpp.clear();
    vaccination_strain_indices_tmp_cpp.clear();

  //  Rcpp::Rcout << "Sec 4: " << i << std::endl;

    infection_history = infection_history_mat(i-1, _);
  //  Rcpp::Rcout << "Sec 4i: " << i << std::endl;

    //  vaccination_history = vaccination_history_mat(i-1,_);

    indices = infection_history > 0;
   // Rcpp::Rcout << "Sec 4ii: " << i << std::endl;

    infection_times = circulation_times[indices];

   // Rcpp::Rcout << "Sec 4iii: " << i << std::endl;

  if (vac_flag) {
      for (int k = 0; k < Ni; k++) { 
        indiv_indic = individuals_vacc_vec[k] == i; // Get boolean values for the individual
        // Extract the individuals details
        if (indiv_indic == true) {
          vac_virus_vec_ind.push_back(vac_virus_vec[k]);
          vac_time_vec_ind.push_back(vac_time_vec[k]);
          vac_flag_ind_cpp.push_back(vac_flag_vec[k]); // Boolean value for when person received previous vaccination
          prev_vac_ind.push_back(prev_vac_vec[k]);
        }
      }
      // Extract data on previous vacciantions
      for (int k = 0; k < vac_flag_ind_cpp.size(); k++) {
        if (vac_flag_ind_cpp[k] == true) {
         // vaccination_strains.push_back(vac_virus_vec_ind[k]);
          vaccination_times_cpp.push_back(vac_time_vec_ind[k]);
          vaccinations_previous_cpp.push_back(prev_vac_ind[k]);
          vaccination_strain_indices_tmp_cpp.push_back(circulation_times_indices[k]);
        }
      }

      vac_flag_ind = vac_flag_ind_cpp;
      vaccination_times = vaccination_times_cpp;
      vaccinations_previous = vaccinations_previous_cpp;
      vaccination_strain_indices_tmp = vaccination_strain_indices_tmp_cpp;
    //  Rcpp::Rcout << "Sec 4v: " << i << std::endl;

     // Rcpp::Rcout << "Redefine everytime: " << std::endl;

    //  indices_vac_individ = individuals_vacc_vec == i;

      //Rcpp::Rcout << "Sec 4v, Post definitions: " << std::endl;

   //   vac_virus_vec_ind = vac_virus_vec[indices_vac_individ]; // length of total ind
    //  vac_time_vec_ind = vac_time_vec[indices_vac_individ];    // length of total ind
     // vac_flag_ind = vac_flag_vec[indices_vac_individ];     // length of total ind
     // prev_vac_ind = prev_vac_vec[indices_vac_individ];     // length of total ind
      
    //  vac_flag_ind_bool = vac_flag_ind > 0;

     // vaccination_times = vac_time_vec_ind[vac_flag_ind_bool];     // just length of vac
     // vaccinations_previous = prev_vac_ind[vac_flag_ind_bool];
      //indices_vac = vac_flag_ind > 0;
     // Rcpp::Rcout << "Sec 4v, Post definitions: " << std::endl;

   //   vaccination_strain_indices_tmp = circulation_times_indices[vac_flag_ind_bool];
    } else {
      vac_flag_ind_bool = false; 
      vac_flag_ind = 0;
      vaccination_times = 0;
      vaccination_strain_indices_tmp = 0;
      vaccinations_previous = 0;
    }
  
   // Rcpp::Rcout << "Sec 5: " << i << std::endl;
    // Only solve is this individual has had infections or vaccinations
    if (infection_times.size() > 0 || vaccination_times.size() > 0) {
   //   Rcpp::Rcout << "Sec 5: " << i << ". Get the indicies." << std::endl;
      infection_strain_indices_tmp = circulation_times_indices[indices];

      index_in_samples = rows_per_indiv_in_samples[i-1]; // count number of samples taken per individual including runs
      end_index_in_samples = rows_per_indiv_in_samples[i] - 1;
      start_index_in_data = cum_nrows_per_individual_in_data[i-1];

      // ====================================================== //
      // =============== CHOOSE MODEL TO SOLVE =============== //
      // ====================================================== //
      // Go to sub function - this is where we have options for different models
      // Note, these are in "boosting_functions.cpp"
      if (base_function) {
     // Rcpp::Rcout << "Sec 5: " << i << ". Call titre_data_fast_individual_base." << std::endl;
	titre_data_fast_individual_base(predicted_titres, mu, mu_short,
					wane, tau,
                    vac_flag,
               //     vac_flag_ind,
                    mu_vac, mu_short_vac, wane_vac, tau_prev_vac,
					infection_times,
					infection_strain_indices_tmp,
                    vaccination_times,
                    vaccination_strain_indices_tmp,
                    vaccinations_previous,
					measurement_strain_indices,
					sample_times,
					index_in_samples,
					end_index_in_samples,
					start_index_in_data,
					nrows_per_blood_sample,
					number_strains,
					antigenic_map_short,
					antigenic_map_long,
          antigenic_map_short_vac,
					antigenic_map_long_vac,
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
                    vac_flag,
           //         vac_flag_ind,
                   // vac_flag_ind,
                    mu_vac, mu_short_vac, wane_vac, tau_prev_vac,
					infection_times,
					infection_strain_indices_tmp,
                    vaccination_times,
                    vaccination_strain_indices_tmp,
                    vaccinations_previous,
					measurement_strain_indices,
					sample_times,
					index_in_samples,
					end_index_in_samples,
					start_index_in_data,
					nrows_per_blood_sample,
					number_strains,
					antigenic_map_short,
					antigenic_map_long,
          antigenic_map_short_vac,
					antigenic_map_long_vac,
					boost_before_infection);
      }
     
    }
  }
 // Rcpp::Rcout << "End of titre_data_fast (infection_model_fast.cpp)" << std::endl;

  return(predicted_titres);
}

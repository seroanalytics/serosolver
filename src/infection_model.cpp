#include <cmath>
#include "wane_function.h"
#include "boosting_functions.h"
#include "helpers.h"
#include "likelihood_funcs.h"
#include "infection_model.h"

//' Overall model function, fast implementation
//'
//' See documentation for \code{\link{titre_data_group}}, as the interface is almost identical
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
			      Nullable<List> additional_arguments
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

  // Parameters for solving model (waning and seniority)
  double sampling_time;
  double time;
  double n_inf;
  
  // Only use the infections that actually happened
  IntegerVector infection_history(number_strains);
  LogicalVector indices;
  
  NumericVector infection_times;
  IntegerVector infection_strain_indices_tmp;
  
  // Pull out model parameters so only need to allocate once
  double mu = theta["mu"];
  double mu_short = theta["mu_short"];
  double wane = theta["wane"];
  double tau = theta["tau"];
  double wane_amount;
  double seniority;
  
  // To store calculated titres
  NumericVector predicted_titres(total_titres);
  
  // For each individual
  for(int i = 1; i <= n; ++i){
    infection_history = infection_history_mat(i-1,_);
    indices = infection_history > 0;
    infection_times = circulation_times[indices];
    
    // Only solve is this individual has had infections
    if(infection_times.size() > 0){
      infection_strain_indices_tmp = circulation_times_indices[indices];
    
      index_in_samples = rows_per_indiv_in_samples[i-1];
      end_index_in_samples = rows_per_indiv_in_samples[i] - 1;
      number_samples = end_index_in_samples - index_in_samples;      
      start_index_in_data = cum_nrows_per_individual_in_data[i-1];
      
      // Go to sub function - this is where we'd have options for different models
      // Note, this is in "boosting_functions.cpp"
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
  return(predicted_titres);
}


  
// Does some precomputation to speed up the model solving code
// Sets up the waning rate and masked infectoin history vectors
void setup_waning_and_masked_cumulative(const NumericVector &theta,
					const IntegerVector &infection_history,
					IntegerVector &cumu_infection_history, 
					IntegerVector &masked_infection_history,
					NumericVector &waning,
					NumericVector &seniority,
					const NumericVector &infection_times,
					const int &max_infections, 
					const double &sampling_time,
					const int &wane_type,
					const double DOB=0,
					const Nullable<List> &additional_arguments=R_NilValue){
  double circulation_time=infection_times[0];
  double time = sampling_time - circulation_time;
  double wane = theta["wane"];
  double tau = theta["tau"];
  double val;  
  
  /* =========================================
     ======== SET UP WANING RATE VECTOR =========
     ========================================= */ 
  // If not linear
  if(wane_type == 0){
    waning[0] = MAX(0, 1.0-wane*time);// Else if linear
  } else {
    // Calculate the waning at the time since infection
    val= wane_function(theta, time, wane);
    waning[0] = MAX(0, 1.0-val);
  }

  // Seniority
  seniority[0] = 1;

  /* =========================================
     ====== SET UP CUMULATIVE INFECTIONS =====
     ========================================= */   
  if(circulation_time > sampling_time) masked_infection_history[0]=0;
  else masked_infection_history[0] = infection_history[0];
  
  cumu_infection_history[0] = masked_infection_history[0];
  
  for(int i = 1; i < max_infections; ++i){
    circulation_time = infection_times[i];
    time = (sampling_time - circulation_time);

    /* =========================================
       ====== SET UP CUMULATIVE INFECTIONS =====
       ========================================= */
    /* Check if isolation time is after the sampling time.
       if so, then we do not test against this strain */ 
    if(circulation_time > sampling_time){
      masked_infection_history[i] = 0;
    } else {
      masked_infection_history[i] = infection_history[i];
    }
    cumu_infection_history[i] = masked_infection_history[i] + cumu_infection_history[i-1];
    
    /* =========================================
       ======== SET UP WANING RATE VECTOR =========
       ========================================= */ 
    // If not linear
    if(wane_type == 0){
      waning[i] = MAX(0, 1.0-wane*time); // Else linear
    } else {
      val= wane_function(theta, time, wane);
      waning[i] = MAX(0, 1.0-val);
    }

    // Seniority
    seniority[i] = MAX(0, 1.0 - tau*(cumu_infection_history[i]-1));
  }  
}


//' Model function sample
//'
//' The main model solving function for a single individual for a single blood sample.
//' NOTES:
//' - Do we want infection history to be a vector of infection times?
//' - Treat the contents of infection_history as a parameter (ie. exposure type)
//' @param theta NumericVector, the named vector of model parameters
//' @param infection_history IntegerVector, the vector of 1s and 0s showing presence/absence of infection for each possible time. 
//' @param infection_times NumericVector, the actual times of circulation that the infection history vector corresponds to
//' @param infection_times_indices IntegerVector, which entry in the melted antigenic map that these infection times correspond to
//' @param samplingTime double, the real time that the sample was taken
//' @param measurement_strain_indices IntegerVector, the indices of all measured strains in the melted antigenic map
//' @param antigenic_map_long NumericVector, the collapsed cross reactivity map for long term boosting, after multiplying by sigma1
//' @param antigenic_map_short NumericVector, the collapsed cross reactivity map for short term boosting, after multiplying by sigma2
//' @param number_strains int, the maximum number of infections that an individual could experience
//' @param DOB double, the date of birth of this individual. Currently not used.
//' @param additional_arguments, Nullable<List> the idea is to use this object to pass more flexible additional arguments to the bottom of the call stack (used for strain dependent boosting right now)
//' @return NumericVector of predicted titres for each entry in measurement_strain_indices
//' @useDynLib serosolver
//' @export
//' @family titre_model
//[[Rcpp::export]]
NumericVector infection_model_indiv(const NumericVector &theta, // Parameter vector
				    const IntegerVector &infection_history, // vector of 1s and 0s for infections
				    const NumericVector &infection_times, // Time of these infections
				    const IntegerVector &infection_times_indices, // Where these infection times fit in antigenic map
				    const double &samplingTime,  // This sampling time
				    const IntegerVector &measurement_strain_indices, // Indices of measured strains in antigenic map
				    const NumericVector &antigenic_map_long,
				    const NumericVector &antigenic_map_short, 
				    const int &number_strains, // Maximum number of infections an individual could experience, if alive the whole time
				    const double &DOB=0,
				    const Nullable<List> &additional_arguments=R_NilValue
				    ){
  double circulation_time;
  
  // Which waning function type should we use
  // 0 is linear decrease
  // 1 is piecewise linear
  int wane_type = theta["wane_type"]; 
  // If titre dependent boosting or not
  bool titre_dependent_boosting = theta["titre_dependent"] == 1;

  // We will need to loop over each strain that was tested
  int n_samples = measurement_strain_indices.size(); // Number of time points sampled
  int max_infections = infection_times.size(); // max number of infections is one for each strain
  
  // Only recording titres for which we have data
  NumericVector predicted_titre(n_samples);
  NumericVector monitored_titres(max_infections);
  // But need to record infection info for each strain that individual could have been
  // infected with
  IntegerVector cum_infection_history(max_infections);
  //IntegerVector cum_infection_history = cumsum(infection_history);

  // Trying just ignoring this, as subset in function above
  IntegerVector masked_infection_history(max_infections); // To avoid solving the model for infections that didn't happen
  NumericVector waning(max_infections);
  NumericVector seniority(max_infections);
  
   // Set up cumulative infection history, masked infection history and waning vector
  setup_waning_and_masked_cumulative(theta,
				     infection_history,
				     cum_infection_history, 
				     masked_infection_history, 
				     waning,	
				     seniority,			     
				     infection_times, 
				     max_infections, 
				     samplingTime,
				     wane_type, 
				     DOB,
				     additional_arguments);

  add_multiple_infections_boost(predicted_titre, monitored_titres,
				theta, 
				infection_times, cum_infection_history, 
				masked_infection_history, 
				infection_times_indices, measurement_strain_indices, 
				antigenic_map_long, antigenic_map_short,
				waning, 
				seniority,
				number_strains, 
				n_samples,
				max_infections,
				titre_dependent_boosting, 
				DOB,
				additional_arguments
				);
  return(predicted_titre);
}

//' Model function individual
//'
//' The main model solving function for a single individual for a vector of sampling times
//' @param theta NumericVector, the named vector of model parameters
//' @param infection_history IntegerVector, the vector of 1s and 0s showing presence/absence of infection for each possible time. 
//' @param circulation_times NumericVector, the actual times of circulation that the infection history vector corresponds to
//' @param circulationMapIndices IntegerVector, which entry in the melted antigenic map that these infection times correspond to
//' @param sample_times NumericVector, the times that each blood sample was taken
//' @param dataIndices IntegerVector, the indices in the overall titre data vector (of observations) that each sample corresponds to
//' @param measurement_strain_indices IntegerVector, the indices of all measured strains in the melted antigenic map
//' @param antigenic_map_long NumericVector, the collapsed cross reactivity map for long term boosting, after multiplying by sigma1
//' @param antigenic_map_short NumericVector, the collapsed cross reactivity map for short term boosting, after multiplying by sigma2
//' @param number_strains int, the maximum number of infections that an individual could experience
//' @param DOB double, the date of birth of this individual. Currently not used.
//' @param additional_arguments, Nullable<List> currently not used, but the idea is to use thsi object to pass more flexible additional arguments to the bottom of the call stack
//' @return NumericVector of predicted titres for each entry in measurement_strain_indices
//' @export
//' @family titre_model
//[[Rcpp::export]]
NumericVector titre_data_individual(const NumericVector &theta, 
				    const IntegerVector &infection_history, 
				    const NumericVector &circulation_times, 
				    const IntegerVector &circulation_times_indices,
				    const NumericVector &sample_times,
				    const IntegerVector &data_indices,
				    const IntegerVector &measurement_strain_indices, 
				    const NumericVector &antigenic_map_long, 
				    const NumericVector &antigenic_map_short,
				    const int &number_strains,
				    const double &DOB=0,
				    const Nullable<List>& additional_arguments=R_NilValue){
  int number_samples = sample_times.size();
  int number_measured_strains = measurement_strain_indices.size();
  NumericVector titres(number_measured_strains);
  
  int start_index = 0;
  int end_index = 0;

  LogicalVector indices = infection_history > 0;

  IntegerVector concise_inf_hist = infection_history[indices];
  NumericVector infection_times = circulation_times[indices];
  IntegerVector inf_map_indices = circulation_times_indices[indices];

  IntegerVector tmp_range;

  for(int i = 0; i < number_samples; ++i){
    end_index = start_index + data_indices[i] - 1;
    // Range index twice
    tmp_range = Range(start_index, end_index);
    titres[tmp_range] = infection_model_indiv(theta,concise_inf_hist,infection_times,inf_map_indices,
					       sample_times[i],measurement_strain_indices[tmp_range],
					       antigenic_map_long, antigenic_map_short,number_strains,
					       DOB, additional_arguments);
    start_index = end_index + 1;
  }
  return(titres);
}

//' Model function overall
//'
//' The main model solving function for a single individual for a vector of sampling times
//' @param theta NumericVector, the named vector of model parameters
//' @param infection_history_mat IntegerMatrix, the matrix of 1s and 0s showing presence/absence of infection for each possible time for each individual. 
//' @param circulation_times NumericVector, the actual times of circulation that the infection history vector corresponds to
//' @param circulation_times_indices IntegerVector, which entry in the melted antigenic map that these infection times correspond to
//' @param sample_times NumericVector, the times that each blood sample was taken
//' @param rows_per_indiv_in_samples IntegerVector, Split the sample times and runs for each individual
//' @param cum_nrows_per_individual_in_data IntegerVector, How many rows in the titre data correspond to each individual?
//' @param rows_per_indiv_in_samples IntegerVector, How many rows in titre data correspond to each individual, sample and repeat?
//' @param measurement_strain_indices IntegerVector, the indices of all measured strains in the melted antigenic map
//' @param antigenic_map_long NumericVector, the collapsed cross reactivity map for long term boosting, after multiplying by sigma1 see \code{\link{create_cross_reactivity_vector}}
//' @param antigenic_map_short NumericVector, the collapsed cross reactivity map for short term boosting, after multiplying by sigma2, see \code{\link{create_cross_reactivity_vector}}
//' @param DOBs NumericVector, the date of birth of all individuals. Currently not used.
//' @param additional_arguments, Nullable<List> the idea is to use thsi object to pass more flexible additional arguments to the bottom of the call stack
//' @return NumericVector of predicted titres for each entry in measurement_strain_indices
//' @export
//' @family titre_model
//[[Rcpp::export]]
NumericVector titre_data_group(const NumericVector &theta, 
			       const IntegerMatrix &infection_history_mat, 
			       const NumericVector &circulation_times,
			       const IntegerVector &circulation_times_indices,
			       const NumericVector &sample_times,
			       const IntegerVector &rows_per_indiv_in_samples,  
			       const IntegerVector &cum_nrows_per_individual_in_data,  
			       const IntegerVector &nrows_per_blood_sample,  
			       const IntegerVector &measurement_strain_indices, 
			       const NumericVector &antigenic_map_long, 
			       const NumericVector &antigenic_map_short,
			       const NumericVector &DOBs,
			       const Nullable<List> &additional_arguments=R_NilValue
			       ){
  int n = infection_history_mat.nrow();
  int n_strains = infection_history_mat.ncol();

  NumericVector titres(measurement_strain_indices.size());
  
  bool check_additional_arguments = additional_arguments.isNotNull(); // Precompute this check so only have to do it once
  bool titre_dependent_boosting = theta["titre_dependent"] == 1;

  // These indices allow us to step through the titre data vector
  // as if it were a matrix ie. number of rows for each individual
  // at a time
  int start_index_samples;
  int end_index_samples;
  int start_index_data;
  int end_index_data;
  IntegerVector tmp_range;
  IntegerVector tmp_range_samples;
  int DOB;

  for(int i=1; i <= n; ++i){
    start_index_samples = rows_per_indiv_in_samples[i-1];
    end_index_samples = rows_per_indiv_in_samples[i] - 1;

    start_index_data = cum_nrows_per_individual_in_data[i-1];
    end_index_data = cum_nrows_per_individual_in_data[i] - 1;

    DOB = DOBs[i-1];

    tmp_range = Range(start_index_data, end_index_data);
    tmp_range_samples = Range(start_index_samples, end_index_samples);

    titres[tmp_range] = titre_data_individual(theta,   // Vector of named model parameters
					      infection_history_mat(i-1,_), // Vector of infection history for individual i
					      circulation_times, // Vector of all virus circulation times, same length and ncol infection_history_mat
					      circulation_times_indices, // Gives the corresponding index in the antigenic map vector
					      sample_times[tmp_range_samples],  // Get sampling times for this individual
					      nrows_per_blood_sample[tmp_range_samples], // The 
					      measurement_strain_indices[tmp_range],  // For the indices in the antigenic map to which each titre corresponds
					      antigenic_map_long, 
					      antigenic_map_short,  
					      n_strains,
					      DOB, 
					      additional_arguments); // The total number of strains that circulated
    
  }
  return(titres);
}


//' Likelihood for one individual
//' 
//' See \code{\link{infection_model_indiv}}, does the same thing but returns the log likelihood given some data. This is used by infection_history_proposal_gibbs
//' @export
//[[Rcpp::export]]
double likelihood_data_individual(const NumericVector &theta, 
				  const IntegerVector &infection_history, 
				  const NumericVector &circulation_times, 
				  const IntegerVector &circulation_times_indices,
				  const NumericVector &sample_times,
				  const IntegerVector &data_indices,
				  const IntegerVector &measurement_strain_indices, 
				  const NumericVector &antigenic_map_long, 
				  const NumericVector &antigenic_map_short,
				  const int &number_strains,
				  const NumericVector &data,
				  const NumericVector &titre_shifts,
				  const double &DOB,
				  const Nullable<List> &additional_arguments
				  ){
  int number_samples = sample_times.size();
  int number_measured_strains = measurement_strain_indices.size();
  double lnlike=0;
  NumericVector titres(number_measured_strains);

  int start_index = 0;
  int end_index = 0;

  LogicalVector indices = infection_history > 0;

  IntegerVector concise_inf_hist = infection_history[indices];
  NumericVector infection_times = circulation_times[indices];
  IntegerVector inf_map_indices = circulation_times_indices[indices];

  for(int i = 0; i < number_samples; ++i){ 
    end_index = start_index + data_indices[i] - 1;
    titres[Range(start_index, end_index)] = infection_model_indiv(theta,
								concise_inf_hist,
								infection_times,
								inf_map_indices,
								sample_times[i],
								measurement_strain_indices[Range(start_index,end_index)],
								antigenic_map_long, antigenic_map_short,
								number_strains, DOB, additional_arguments);
    start_index = end_index + 1;
  }
  lnlike = likelihood_titre_basic(titres, data, theta, titre_shifts);
  return(lnlike);
}

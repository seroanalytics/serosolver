#include "antibody_models_individual.h"

#ifndef MAX
#define MAX(a,b) ((a) < (b) ? (b) : (a)) // define MAX function for use later
#endif

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b)) // define MAX function for use later
#endif

void antibody_data_model_individual_new(NumericVector &predicted_antibody_levels,
                                         const NumericVector &start_antibody_levels,
                                         const NumericVector &births,
                                         const double &boost_long,
                                         const double &boost_short,
                                         const double &boost_delay,
                                         const double &wane_short,
                                         const double &wane_long,
                                         const double &antigenic_seniority,
                                         const NumericVector &possible_exposure_times,
                                         const IntegerVector &exposure_indices,
                                         const IntegerVector &biomarker_id_indices,
                                         const IntegerVector &start_level_indices,
                                         const NumericVector &sample_times,
                                         const int &index_in_samples,
                                         const int &end_index_in_samples,
                                         const int &start_index_in_data1,
                                         const IntegerVector &nrows_per_blood_sample,
                                         const int &number_possible_exposures,
                                         const double *antigenic_map_short,
                                         const double *antigenic_map_long,
                                         bool boost_before_infection = false
 ){
   double sampling_time;
   double time;
   double n_inf;
   double wane_short_amount;
   double wane_long_amount;
   double seniority;
   int n_measurements;
   int max_infections = possible_exposure_times.size();
   int end_index_in_data;
   int tmp_measurement_index;
   int start_index_in_data = start_index_in_data1;
   int inf_map_index;
   int index;
   
   // For each sample this individual has
   for(int j = index_in_samples; j <= end_index_in_samples; ++j){
     sampling_time = sample_times[j];
     n_inf = 1.0;
     // Find number of measurements in the predicted_antibody_levels vector that correspond to this sample
     n_measurements = nrows_per_blood_sample[j];
     // Only iterate through indices for this sample
     end_index_in_data = start_index_in_data + n_measurements;
     tmp_measurement_index = start_index_in_data;
     
     // Include starting titre contributions
     // Time elapsed since birth/first time
     // Assume that first entry for birth does not change
     time = sampling_time - births[start_index_in_data];
     wane_long_amount= MAX(0, 1.0 - (wane_long*time)); 
     

     // For each measured marker, find the biomarker id index which will match an entry in start_antibody_levels
     // Add this to the predicted antibody level, with waning
     for(int k = 0; k < n_measurements; ++k){
       index = start_level_indices[tmp_measurement_index + k];
       predicted_antibody_levels[tmp_measurement_index + k] += start_antibody_levels[index]*wane_long_amount;
     }
     
     // Sum all infections that would contribute towards observed antibody levels at this time
     for(int x = 0; x < max_infections; ++x){
       // Only go further if this sample happened after the infection
       if((boost_before_infection && sampling_time > (possible_exposure_times[x] + boost_delay)) ||
          (!boost_before_infection && sampling_time >= (possible_exposure_times[x] + boost_delay))){
         time = sampling_time - (possible_exposure_times[x] + boost_delay); // Time between sample and infection + boost
         wane_short_amount= MAX(0, 1.0 - (wane_short*time)); // Waning of the short-term response
         wane_long_amount= MAX(0, 1.0 - (wane_long*time)); // Waning of the long-term response
         
         seniority = MAX(0, 1.0 - antigenic_seniority*(n_inf - 1.0)); // Antigenic seniority
         inf_map_index = exposure_indices[x]; // Index of this infecting antigen in antigenic map
         
         // Find contribution to each measured antibody level from this infection
         for(int k = 0; k < n_measurements; ++k){
           index = biomarker_id_indices[tmp_measurement_index + k]*number_possible_exposures + inf_map_index;
           predicted_antibody_levels[tmp_measurement_index + k] += seniority *
             ((boost_long*antigenic_map_long[index])*wane_long_amount + (boost_short*antigenic_map_short[index])*wane_short_amount);
         }
         
         ++n_inf;
       }
     }
     start_index_in_data = end_index_in_data;
   }
 }
 
 // [[Rcpp::export]]
Rcpp::NumericVector antibody_model_individual_wrapper(const double &boost_long,
                                                      const double &boost_short,
                                                      const double &boost_delay,
                                                      const double &wane_short,
                                                      const double &wane_long,
                                                      const double &antigenic_seniority,
                                                      const int &birth,
                                                      const NumericVector &start_antibody_levels,
                                                      const int number_possible_exposures,
                                                      const NumericVector &possible_exposure_times,
                                                      const IntegerVector &exposure_indices,
                                                      const IntegerVector &biomarker_id_indices,
                                                      const NumericVector &sample_times,
                                                      const arma::mat &antigenic_map_long,
                                                      const arma::mat &antigenic_map_short){
  
  Rcpp::NumericVector predicted_antibody_levels(biomarker_id_indices.size()*sample_times.size());
  int index_in_samples = 0;
  int end_index_in_samples = index_in_samples + sample_times.size() - 1;
  int start_index_in_data = 0;
  Rcpp::IntegerVector nrows_per_blood_sample(sample_times.size());
  for(int i = 0; i < nrows_per_blood_sample.size(); ++i){
    nrows_per_blood_sample[i] = biomarker_id_indices.size();
  }
  IntegerVector biomarker_id_indices_use = rep(biomarker_id_indices, sample_times.size());
  IntegerVector start_level_indices = biomarker_id_indices_use;
    antibody_data_model_individual_new(predicted_antibody_levels, start_antibody_levels,birth,
                                     boost_long, boost_short, boost_delay, wane_short,wane_long, antigenic_seniority,
                                 possible_exposure_times,exposure_indices, biomarker_id_indices_use, start_level_indices,
                                 sample_times,
                                 index_in_samples, end_index_in_samples, start_index_in_data, nrows_per_blood_sample,number_possible_exposures,
                                 antigenic_map_short.colptr(0), antigenic_map_long.colptr(0),
                                 false);
  
  return(predicted_antibody_levels);
}


//' Antibody model for one individual
//' 
//' A fast implementation of the basic boosting function, giving predicted antibody_levels for a number of samples for one individual. Note that this version attempts to minimise memory allocations.
//' @family antibody_models
//' @export
//' @seealso \code{\link{antibody_model}}
void antibody_data_model_individual(
    NumericVector &predicted_antibody_levels,
				     const double &boost_long,
				     const double &boost_short,
				     const double &boost_delay,
				     const double &wane_short,
				     const double &wane_long,
				     const double &antigenic_seniority,
				     const NumericVector &possible_exposure_times,
				     const IntegerVector &exposure_indices,
				     const IntegerVector &biomarker_id_indices,
				     const NumericVector &sample_times,
				     const int &index_in_samples,
				     const int &end_index_in_samples,
				     const int &start_index_in_data1,
				     const IntegerVector &nrows_per_blood_sample,
				     const int &number_possible_exposures,
				     const double *antigenic_map_short,
				     const double *antigenic_map_long,
				     bool boost_before_infection = false
				     ){
  double sampling_time;
  double time;
  double n_inf;
  double wane_short_amount;
  double wane_long_amount;
  double seniority;
  int n_measurements;
  int max_infections = possible_exposure_times.size();
  int end_index_in_data;
  int tmp_measurement_index;
  int start_index_in_data = start_index_in_data1;
  int inf_map_index;
  int index;

  // For each sample this individual has
  for(int j = index_in_samples; j <= end_index_in_samples; ++j){
    sampling_time = sample_times[j];
    n_inf = 1.0;
    // Find number of measurements in the predicted_antibody_levels vector that correspond to this sample
    n_measurements = nrows_per_blood_sample[j];
    // Only iterate through indices for this sample
    end_index_in_data = start_index_in_data + n_measurements;
    tmp_measurement_index = start_index_in_data;
    // Sum all infections that would contribute towards observed antibody levels at this time
    for(int x = 0; x < max_infections; ++x){
      // Only go further if this sample happened after the infection
      if((boost_before_infection && sampling_time > (possible_exposure_times[x] + boost_delay)) ||
	 (!boost_before_infection && sampling_time >= (possible_exposure_times[x] + boost_delay))){
	time = sampling_time - (possible_exposure_times[x] + boost_delay); // Time between sample and infection + boost
	wane_short_amount= MAX(0, 1.0 - (wane_short*time)); // Waning of the short-term response
	wane_long_amount= MAX(0, 1.0 - (wane_long*time)); // Waning of the long-term response
	
	seniority = MAX(0, 1.0 - antigenic_seniority*(n_inf - 1.0)); // Antigenic seniority
	inf_map_index = exposure_indices[x]; // Index of this infecting antigen in antigenic map

	// Find contribution to each measured antibody level from this infection
for(int k = 0; k < n_measurements; ++k){
	  index = biomarker_id_indices[tmp_measurement_index + k]*number_possible_exposures + inf_map_index;
	  predicted_antibody_levels[tmp_measurement_index + k] += seniority *
	    ((boost_long*antigenic_map_long[index])*wane_long_amount + (boost_short*antigenic_map_short[index])*wane_short_amount);
	}
 
	++n_inf;
      }
    }
    start_index_in_data = end_index_in_data;
  }
}

//' Antibody dependent boosting model, one individual
//' 
//' A fast implementation of the antibody dependent boosting function, giving predicted antibody levels for a number of samples for one individual. Note that this version attempts to minimise memory allocations.
//' @family antibody_models
//' @seealso \code{\link{antibody_model}}
void antibody_dependent_boosting_model_individual(NumericVector &predicted_antibody_levels,
					 const double &boost_long,
					 const double &boost_short,
					 const double &wane_short,
					 const double &antigenic_seniority,
					 const double &gradient,
					 const double &boost_limit,
					 const NumericVector &possible_exposure_times,
					 const IntegerVector &exposure_indices,
					 const IntegerVector &biomarker_id_indices,
					 const NumericVector &sample_times,
					 const int &index_in_samples,
					 const int &end_index_in_samples,
					 const int &start_index_in_data1,
					 const IntegerVector &nrows_per_blood_sample,
					 const int &number_possible_exposures,
					 const double *antigenic_map_short,
					 const double *antigenic_map_long,
				     bool boost_before_infection = false
					 ){
  double sampling_time;
  double time;
  double n_inf;
  double wane_short_amount;
  double seniority;
  double infection_time;

  double boost = 0;
  double long_boost=0;
  double monitored_antibody_level=0;
  double short_boost=0;
  double antibody_level_suppression = MAX(0,1.0 - gradient*boost_limit);

  int n_measurements;
  int max_infections = possible_exposure_times.size();
  int end_index_in_data;
  int tmp_measurement_index;
  int start_index_in_data = start_index_in_data1;
  int inf_map_index;
  int inf_map_index_tmp;
  int index;

  NumericVector monitored_antibody_levels(max_infections);

  // For each sample this individual has
  for(int j = index_in_samples; j <= end_index_in_samples; ++j){
    sampling_time = sample_times[j];
    n_inf = 0.0;

    // Find number of antibody_levels in the predicted_antibody_levels vector that correspond to this sample
    n_measurements = nrows_per_blood_sample[j];

    // Only iterate through indices for this sample
    end_index_in_data = start_index_in_data + n_measurements;
    tmp_measurement_index = start_index_in_data;

    // Sum all infections that would contribute towards observed antibody_levels at this time
    for(int x = 0; x < max_infections; ++x){
      // Only go further if this sample happened after the infection
        if((boost_before_infection && sampling_time > possible_exposure_times[x]) ||
	   (!boost_before_infection && sampling_time >= possible_exposure_times[x])){
	  //      if(sampling_time >= possible_exposure_times[x]){
	  monitored_antibody_level = 0;
	infection_time = possible_exposure_times[x];
	time = sampling_time - infection_time; // Time between sample and infection
	inf_map_index = exposure_indices[x]; // Index of this infecting antigen in antigenic map

	// Add up contribution of all previous infections to antibodies that
	// would be observed at this infection time
	for(int ii = x - 1; ii >= 0; --ii){
	  inf_map_index_tmp = inf_map_index * number_possible_exposures + exposure_indices[ii];
	  seniority = MAX(0, 1.0 - antigenic_seniority * ii);
	  wane_short_amount = MAX(0, 1.0 - wane_short * (infection_time - possible_exposure_times[ii]));

	  long_boost = seniority * boost_long * antigenic_map_long[inf_map_index_tmp];
	  short_boost = seniority * boost_short * antigenic_map_short[inf_map_index_tmp];
	  if(monitored_antibody_levels[ii] >= boost_limit){
	    long_boost *= antibody_level_suppression;
	    short_boost *= antibody_level_suppression;
	  } else {
	    long_boost *= MAX(0,1.0 - gradient*monitored_antibody_levels[ii]);
	    short_boost *= MAX(0,1.0 - gradient*monitored_antibody_levels[ii]);	    
	  }
	  long_boost = MAX(0, long_boost);
	  short_boost = MAX(0, short_boost);
	  boost = long_boost + short_boost * wane_short_amount;
	  monitored_antibody_level += boost;
	}
	monitored_antibody_levels[x] = monitored_antibody_level;

	wane_short_amount= MAX(0, 1.0 - (wane_short*time)); // Basic waning function
	seniority = MAX(0, 1.0 - antigenic_seniority*n_inf); // Antigenic seniority
	
	// Find contribution to each measured antigen from this infection
	for(int k = 0; k < n_measurements; ++k){
	  index = biomarker_id_indices[tmp_measurement_index + k]*number_possible_exposures + inf_map_index;

	  long_boost = seniority * boost_long * antigenic_map_long[index];
	  short_boost = seniority * boost_short * antigenic_map_short[index];
	  
 // Antibody dependent boosting - at ceiling
	  if(monitored_antibody_levels[x] >= boost_limit){
	    long_boost *= antibody_level_suppression;
	    short_boost *= antibody_level_suppression;
// Antibody dependent boosting - below ceiling
	  } else {
	    long_boost = long_boost * (1 - gradient * monitored_antibody_levels[x]); 
	    short_boost = short_boost * (1 - gradient * monitored_antibody_levels[x]); // Antibody dependent boosting - below ceiling
	  }
	  long_boost = MAX(0, long_boost);
	  short_boost = MAX(0, short_boost);
	  boost = long_boost + short_boost * wane_short_amount;
	  predicted_antibody_levels[tmp_measurement_index + k] += boost;
	}
	++n_inf;
      }
    }
    start_index_in_data = end_index_in_data;
  }
}
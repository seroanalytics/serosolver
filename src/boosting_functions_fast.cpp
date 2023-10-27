#include "boosting_functions_fast.h"

#ifndef MAX
#define MAX(a,b) ((a) < (b) ? (b) : (a)) // define MAX function for use later
#endif

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b)) // define MAX function for use later
#endif

//' Base boosting fast
//' 
//' A fast implementation of the basic boosting function, giving predicted titres for a number of samples for one individual. Note that this version attempts to minimise memory allocations.
//' @family boosting_functions
//' @seealso \code{\link{titre_data_fast}}
void titre_data_fast_individual_base(NumericVector &predicted_titres,
				     const double &mu,
				     const double &mu_short,
				     const double &wane,
				     const double &tau,
				     const NumericVector &wane_amounts,
				     const NumericVector &seniority_amounts,
				     const NumericVector &infection_times,
				     const IntegerVector &infection_strain_indices_tmp,
				     const IntegerVector &measurement_strain_indices,
				     const NumericVector &sample_times,
				     const int &index_in_samples,
				     const int &end_index_in_samples,
				     const int &start_index_in_data1,
				     const IntegerVector &nrows_per_blood_sample,
				     const int &number_strains,
				     const NumericVector &antigenic_map_short,
				     const NumericVector &antigenic_map_long,
				     bool boost_before_infection = false
				     ){
  double sampling_time;
  double time;
  int n_inf;
  double wane_amount;
  double seniority;

  int n_titres;
  int max_infections = infection_times.size();
  int end_index_in_data;
  int tmp_titre_index;
  int start_index_in_data = start_index_in_data1;
  int inf_map_index;
  int index;

  // For each sample this individual has
  for(int j = index_in_samples; j <= end_index_in_samples; ++j){
    sampling_time = sample_times[j];
    n_inf = 1;

    // Find number of titres in the predicted_titres vector that correspond to this sample
    n_titres = nrows_per_blood_sample[j];
    // Only iterate through indices for this sample
    end_index_in_data = start_index_in_data + n_titres;
    tmp_titre_index = start_index_in_data;
    // Sum all infections that would contribute towards observed titres at this time
    for(int x = 0; x < max_infections; ++x){
      // Only go further if this sample happened after the infection
      //if(sampling_time >= infection_times[x]){
      if((boost_before_infection && sampling_time > infection_times[x]) ||
	 (!boost_before_infection && sampling_time >= infection_times[x])){
    	time = sampling_time - infection_times[x]; // Time between sample and infection
    	wane_amount= wane_amounts[time]; // MAX(0, 1.0 - (wane*time)); //wane_amounts[time]; //  Basic waning function
    	seniority = seniority_amounts[n_inf-1]; // MAX(0, 1.0 - tau*(n_inf - 1.0)); //  //  Antigenic seniority
    	inf_map_index = infection_strain_indices_tmp[x]; // Index of this infecting strain in antigenic map
    
    	// Find contribution to each measured titre from this infection
    	for(int k = 0; k < n_titres; ++k){
    	  index = measurement_strain_indices[tmp_titre_index + k]*number_strains + inf_map_index;
    	  predicted_titres[tmp_titre_index + k] += seniority *
    	    ((mu*antigenic_map_long[index]) + (mu_short*antigenic_map_short[index])*wane_amount);
    	}
	++n_inf;
      }
    }
    start_index_in_data = end_index_in_data;
  }
}


//' Alternative waning fast
//' 
//' A fast implementation of the alternative waning function, giving predicted titres for a number of samples for one individual. Note that this version attempts to minimise memory allocations.
//' @family boosting_functions
//' @seealso \code{\link{titre_data_fast}}
void titre_data_fast_individual_wane2(NumericVector &predicted_titres,
				      const double &mu,
				      const double &mu_short,
				      const double &wane,
				      const double &tau,
				      const double &kappa,
				      const double &t_change,
				      const NumericVector &infection_times,
				      const IntegerVector &infection_strain_indices_tmp,
				      const IntegerVector &measurement_strain_indices,
				      const NumericVector &sample_times,
				      const int &index_in_samples,
				      const int &end_index_in_samples,
				      const int &start_index_in_data1,
				      const IntegerVector &nrows_per_blood_sample,
				      const int &number_strains,
				      const NumericVector &antigenic_map_short,
				      const NumericVector &antigenic_map_long,
				     bool boost_before_infection = false
				      ){
  double sampling_time;
  double time;
  double n_inf;
  double wane_amount;
  double seniority;

  double wane_2 = -kappa*wane;
  double wane_2_val; // Interaction term

  int n_titres;
  int max_infections = infection_times.size();
  int end_index_in_data;
  int tmp_titre_index;
  int start_index_in_data = start_index_in_data1;
  int inf_map_index;
  int index;

  // For each sample this individual has
  for(int j = index_in_samples; j <= end_index_in_samples; ++j){
    sampling_time = sample_times[j];
    n_inf = 1.0;

    // Find number of titres in the predicted_titres vector that correspond to this sample
    n_titres = nrows_per_blood_sample[j];

    // Only iterate through indices for this sample
    end_index_in_data = start_index_in_data + n_titres;
    tmp_titre_index = start_index_in_data;

    // Sum all infections that would contribute towards observed titres at this time
    for(int x = 0; x < max_infections; ++x){
      // Only go further if this sample happened after the infection
      //if(sampling_time >= infection_times[x]){
      if((boost_before_infection && sampling_time > infection_times[x]) ||
	 (!boost_before_infection && sampling_time >= infection_times[x])){
	time = sampling_time - infection_times[x]; // Time between sample and infection

	/////////////////////////////////
	// Amanda's alternative waning function
	if(time > t_change){
	  wane_2_val = wane_2*(time - t_change); 
	}else{
	  wane_2_val = 0;
	}
	wane_amount = MAX(0, 1.0 - (wane*time+wane_2_val));
	/////////////////////////////////

	seniority = MAX(0, 1.0 - tau*(n_inf - 1.0)); // Antigenic seniority
	inf_map_index = infection_strain_indices_tmp[x]; // Index of this infecting strain in antigenic map

	// Find contribution to each measured titre from this infection
	for(int k = 0; k < n_titres; ++k){
	  index = measurement_strain_indices[tmp_titre_index + k]*number_strains + inf_map_index;
	  predicted_titres[tmp_titre_index + k] += seniority *
	    ((mu*antigenic_map_long[index]) + (mu_short*antigenic_map_short[index])*wane_amount);
	}
	++n_inf;
      }
    }
    start_index_in_data = end_index_in_data;
  }
}


//' Titre dependent boosting fast
//' 
//' A fast implementation of the titre dependent boosting function, giving predicted titres for a number of samples for one individual. Note that this version attempts to minimise memory allocations.
//' @family boosting_functions
//' @seealso \code{\link{titre_data_fast}}
void titre_data_fast_individual_titredep(NumericVector &predicted_titres,
					 const double &mu,
					 const double &mu_short,
					 const double &wane,
					 const double &tau,
					 const double &gradient,
					 const double &boost_limit,
					 const NumericVector &infection_times,
					 const IntegerVector &infection_strain_indices_tmp,
					 const IntegerVector &measurement_strain_indices,
					 const NumericVector &sample_times,
					 const int &index_in_samples,
					 const int &end_index_in_samples,
					 const int &start_index_in_data1,
					 const IntegerVector &nrows_per_blood_sample,
					 const int &number_strains,
					 const NumericVector &antigenic_map_short,
					 const NumericVector &antigenic_map_long,
				     bool boost_before_infection = false
					 ){
  double sampling_time;
  double time;
  double n_inf;
  double wane_amount;
  double seniority;
  double infection_time;

  double boost = 0;
  double long_boost=0;
  double monitored_titre=0;
  double short_boost=0;
  double titre_suppression = MAX(0,1.0 - gradient*boost_limit);

  int n_titres;
  int max_infections = infection_times.size();
  int end_index_in_data;
  int tmp_titre_index;
  int start_index_in_data = start_index_in_data1;
  int inf_map_index;
  int inf_map_index_tmp;
  int index;

  NumericVector monitored_titres(max_infections);

  // For each sample this individual has
  for(int j = index_in_samples; j <= end_index_in_samples; ++j){
    sampling_time = sample_times[j];
    n_inf = 0.0;

    // Find number of titres in the predicted_titres vector that correspond to this sample
    n_titres = nrows_per_blood_sample[j];

    // Only iterate through indices for this sample
    end_index_in_data = start_index_in_data + n_titres;
    tmp_titre_index = start_index_in_data;

    // Sum all infections that would contribute towards observed titres at this time
    for(int x = 0; x < max_infections; ++x){
      // Only go further if this sample happened after the infection
        if((boost_before_infection && sampling_time > infection_times[x]) ||
	   (!boost_before_infection && sampling_time >= infection_times[x])){
	  //      if(sampling_time >= infection_times[x]){
	monitored_titre = 0;
	infection_time = infection_times[x];
	time = sampling_time - infection_time; // Time between sample and infection
	inf_map_index = infection_strain_indices_tmp[x]; // Index of this infecting strain in antigenic map

	// Add up contribution of all previous infections to titre that
	// would be observed at this infection time
	for(int ii = x - 1; ii >= 0; --ii){
	  inf_map_index_tmp = inf_map_index * number_strains + infection_strain_indices_tmp[ii];
	  seniority = MAX(0, 1.0 - tau * ii);
	  wane_amount = MAX(0, 1.0 - wane * (infection_time - infection_times[ii]));

	  long_boost = seniority * mu * antigenic_map_long[inf_map_index_tmp];
	  short_boost = seniority * mu_short * antigenic_map_short[inf_map_index_tmp];
	  if(monitored_titres[ii] >= boost_limit){
	    long_boost *= titre_suppression;
	    short_boost *= titre_suppression;
	  } else {
	    long_boost *= MAX(0,1.0 - gradient*monitored_titres[ii]);
	    short_boost *= MAX(0,1.0 - gradient*monitored_titres[ii]);	    
	  }
	  long_boost = MAX(0, long_boost);
	  short_boost = MAX(0, short_boost);
	  boost = long_boost + short_boost * wane_amount;
	  monitored_titre += boost;
	}
	monitored_titres[x] = monitored_titre;

	wane_amount= MAX(0, 1.0 - (wane*time)); // Basic waning function
	seniority = MAX(0, 1.0 - tau*n_inf); // Antigenic seniority
	
	// Find contribution to each measured titre from this infection
	for(int k = 0; k < n_titres; ++k){
	  index = measurement_strain_indices[tmp_titre_index + k]*number_strains + inf_map_index;

	  long_boost = seniority * mu * antigenic_map_long[index];
	  short_boost = seniority * mu_short * antigenic_map_short[index];
	  
 // Titre dependent boosting - at ceiling
	  if(monitored_titres[x] >= boost_limit){
	    long_boost *= titre_suppression;
	    short_boost *= titre_suppression;
// Titre dependent boosting - below ceiling
	  } else {
	    long_boost = long_boost * (1 - gradient * monitored_titres[x]); 
	    short_boost = short_boost * (1 - gradient * monitored_titres[x]); // Titre dependent boosting - below ceiling
	  }
	  long_boost = MAX(0, long_boost);
	  short_boost = MAX(0, short_boost);
	  boost = long_boost + short_boost * wane_amount;
	  predicted_titres[tmp_titre_index + k] += boost;
	}
	++n_inf;
      }
    }
    start_index_in_data = end_index_in_data;
  }
}



//' Base boosting fast
//' 
//' A fast implementation of the basic boosting function, giving predicted titres for a number of samples for one individual. Note that this version attempts to minimise memory allocations.
//' @family boosting_functions
//' @seealso \code{\link{titre_data_fast}}
void titre_data_fast_individual_strain_dependent(NumericVector &predicted_titres,
						 const NumericVector &mus,
						 const IntegerVector &boosting_vec_indices,
						 const double &mu_short,
						 const double &wane,
						 const double &tau,
						 const NumericVector &infection_times,
						 const IntegerVector &infection_strain_indices_tmp,
						 const IntegerVector &measurement_strain_indices,
						 const NumericVector &sample_times,
						 const int &index_in_samples,
						 const int &end_index_in_samples,
						 const int &start_index_in_data1,
						 const IntegerVector &nrows_per_blood_sample,
						 const int &number_strains,
						 const NumericVector &antigenic_map_short,
						 const NumericVector &antigenic_map_long,
				     bool boost_before_infection = false
						 ){
  double sampling_time;
  double time;
  double n_inf;
  double wane_amount;
  double seniority;

  double mu;

  int n_titres;
  int max_infections = infection_times.size();
  int end_index_in_data;
  int tmp_titre_index;
  int start_index_in_data = start_index_in_data1;
  int inf_map_index;
    int index;
  
  // For each sample this individual has
  for(int j = index_in_samples; j <= end_index_in_samples; ++j){
    sampling_time = sample_times[j];
    n_inf = 1.0;

    // Find number of titres in the predicted_titres vector that correspond to this sample
    n_titres = nrows_per_blood_sample[j];

    // Only iterate through indices for this sample
    end_index_in_data = start_index_in_data + n_titres;
    tmp_titre_index = start_index_in_data;

    // Sum all infections that would contribute towards observed titres at this time
    for(int x = 0; x < max_infections; ++x){
      // Only go further if this sample happened after the infection
        if((boost_before_infection && sampling_time > infection_times[x]) ||
	   (!boost_before_infection && sampling_time >= infection_times[x])){
	  //if(sampling_time >= infection_times[x]){
	time = sampling_time - infection_times[x]; // Time between sample and infection
	wane_amount= MAX(0, 1.0 - (wane*time)); // Basic waning function
	seniority = MAX(0, 1.0 - tau*(n_inf - 1.0)); // Antigenic seniority
	inf_map_index = infection_strain_indices_tmp[x]; // Index of this infecting strain in antigenic map
	mu = mus[boosting_vec_indices[inf_map_index]];
	// Find contribution to each measured titre from this infection
	for(int k = 0; k < n_titres; ++k){
	  index = measurement_strain_indices[tmp_titre_index + k]*number_strains + inf_map_index;
	  predicted_titres[tmp_titre_index + k] += seniority *
	    ((mu*antigenic_map_long[index]) + (mu_short*antigenic_map_short[index])*wane_amount);
	}
	++n_inf;
      }
    }
    start_index_in_data = end_index_in_data;
  }
}

//' Base boosting fast
 //' 
 //' A fast implementation of the basic boosting function, giving predicted titres for a number of samples for one individual. Note that this version attempts to minimise memory allocations.
 //' @family boosting_functions
 //' @seealso \code{\link{titre_data_fast}}
 void titre_data_fast_individual_complex_cr(NumericVector &predicted_titres,
                                      const double &mu,
                                      const double &mu_short,
                                      const double &wane,
                                      const double &tau,
                                      
                                      // Cross reactivity parameters
                                      const double &sigma_s, // short term cross reactivity
                                      const double &sigma_l, // long term cross reactivity
                                      const double &sigma_birth_mod_s, // short term cross-reactivity modifier to pre-birth
                                      const double &sigma_birth_mod_l, // long term cross-reactivity modifier to pre-birth
                                      const double &sigma_future_mod_s, // short term cross-reactivity modifier to future
                                      const double &sigma_future_mod_l, // long term cross-reactivity modifier to future
                                      
                                      const NumericVector &wane_amounts,
                                      const NumericVector &seniority_amounts,
                                      const NumericVector &infection_times,
                                      const IntegerVector &infection_strain_indices_tmp,
                                      const IntegerVector &measurement_strain_indices,
                                      const NumericVector &sample_times,
                                      const int &index_in_samples,
                                      const int &end_index_in_samples,
                                      const int &start_index_in_data1,
                                      const IntegerVector &nrows_per_blood_sample,
                                      const int &number_strains,
                                      const NumericVector &antigenic_distances,
                                      bool boost_before_infection = false
 ){
     double sampling_time;
     double time;
     int n_inf;
     double wane_amount;
     double seniority;
     
     double cr_short;
     double cr_long;
     int n_titres;
     int max_infections = infection_times.size();
     int end_index_in_data;
     int tmp_titre_index;
     int start_index_in_data = start_index_in_data1;
     int inf_map_index;

     // Indices in antigenic map of first and last infecting strain
     int first_strain = infection_strain_indices_tmp[0];
     int last_strain;// = infection_strain_indices_tmp[max_infections-1];
     
     double intercept_long_future;
     double intercept_short_future;
     
     double intercept_long_past;
     double intercept_short_past;
     
     int last_strain_index;
     int measured_strain_with_last_index;
     int infecting_strain_with_first_index;
     int infecting_strain_with_last_index;
     int measured_strain_with_first_index;
     int measured_strain;
     int measured_strain_index;

     // For each sample this individual has
     for(int j = index_in_samples; j <= end_index_in_samples; ++j){
         sampling_time = sample_times[j];
         n_inf = 1;
         // Find number of titres in the predicted_titres vector that correspond to this sample
         n_titres = nrows_per_blood_sample[j];
         // Only iterate through indices for this sample
         end_index_in_data = start_index_in_data + n_titres;
         tmp_titre_index = start_index_in_data;
         
         // Find last infection on or before this sampling time
         last_strain = max_infections - 1;
         while(infection_times[last_strain] > sampling_time && last_strain > 0){
             last_strain--;
         }
         last_strain_index = infection_strain_indices_tmp[last_strain];

         // Sum all infections that would contribute towards observed titres at this time
         for(int x = 0; x < max_infections; ++x){
             // Only go further if this sample happened after the infection
             //if(sampling_time >= infection_times[x]){
             if((boost_before_infection && sampling_time > infection_times[x]) ||
                (!boost_before_infection && sampling_time >= infection_times[x])){
                 time = sampling_time - infection_times[x]; // Time between sample and infection
                 wane_amount= wane_amounts[time]; // MAX(0, 1.0 - (wane*time)); // //  Basic waning function
                 seniority = seniority_amounts[n_inf-1];//MAX(0, 1.0 - tau*(n_inf - 1.0)); //  //  Antigenic seniority
                 
                 inf_map_index = infection_strain_indices_tmp[x]; // Index of this infecting strain in antigenic map
                 
                 // Infecting strain with last strain
                 infecting_strain_with_last_index = last_strain_index*number_strains + inf_map_index;
                 // Infecting strain with first strain
                 infecting_strain_with_first_index = inf_map_index*number_strains + first_strain;
                 
                 intercept_long_future = MAX(0, 1.0 - sigma_l*antigenic_distances[infecting_strain_with_last_index]);
                 intercept_short_future = MAX(0, 1.0 - sigma_s*antigenic_distances[infecting_strain_with_last_index]);
                 
                 intercept_long_past = MAX(0, 1.0 - sigma_l*antigenic_distances[infecting_strain_with_first_index]);
                 intercept_short_past = MAX(0, 1.0 - sigma_s*antigenic_distances[infecting_strain_with_first_index]);
                 
                 // Find contribution to each measured titre from this infection
                 for(int k = 0; k < n_titres; ++k){
                     // Last strain with measured strain
                     measured_strain_with_last_index = measurement_strain_indices[tmp_titre_index + k]*number_strains + last_strain_index;
                     // First strain with measured strain
                     measured_strain_with_first_index = measurement_strain_indices[tmp_titre_index + k]*number_strains + first_strain;
                     
                     // Infecting strain with measured strain
                     measured_strain = measurement_strain_indices[tmp_titre_index + k];
                     measured_strain_index = measurement_strain_indices[tmp_titre_index + k]*number_strains + inf_map_index;
                     
                     if(measured_strain > last_strain_index){
                         cr_long = MAX(0, intercept_long_future - (sigma_l + sigma_future_mod_l)*antigenic_distances[measured_strain_with_last_index]);
                         cr_short = MAX(0, intercept_short_future - (sigma_s + sigma_future_mod_s)*antigenic_distances[measured_strain_with_last_index]);
                     
                     } else if(measured_strain < first_strain) {
                         cr_long = MAX(0, intercept_long_past - (sigma_l + sigma_birth_mod_l)*antigenic_distances[measured_strain_with_first_index]);
                         cr_short = MAX(0, intercept_short_past - (sigma_s + sigma_birth_mod_s)*antigenic_distances[measured_strain_with_first_index]);
                     } else {
                         cr_long = MAX(0, 1.0 - (sigma_l)*antigenic_distances[measured_strain_index]);
                         cr_short = MAX(0, 1.0 - (sigma_s)*antigenic_distances[measured_strain_index]);
                     }
                     predicted_titres[tmp_titre_index + k] += seniority *
                         ((mu*cr_long) + mu_short*cr_short*wane_amount);
                 }
                 ++n_inf;
             }
         }
         start_index_in_data = end_index_in_data;
     }
 }

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
	wane_amount= MAX(0, 1.0 - (wane*time)); // Basic waning function
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

  //NumericVector monitored_titres(max_infections);
 // Rcpp::Rcout << std::endl << std::endl << "Starting" << std::endl;
  // For each sample this individual has
  for(int j = index_in_samples; j <= end_index_in_samples; ++j){
    sampling_time = sample_times[j];
   // Rcpp::Rcout << "Sampling time: " << sampling_time << std::endl;
    n_inf = 0.0;

    // Find number of titres in the predicted_titres vector that correspond to this sample
    n_titres = nrows_per_blood_sample[j];

    // For each measured strain, find the titre at the time of each infection 
    NumericMatrix monitored_titres(max_infections, n_titres);
    
    // Only iterate through indices for this sample
    end_index_in_data = start_index_in_data + n_titres;
    tmp_titre_index = start_index_in_data;

    // Sum all infections that would contribute towards observed titres at this time
    for(int x = 0; x < max_infections; ++x){
      // Only go further if this sample happened after the infection
        if((boost_before_infection && sampling_time > infection_times[x]) ||
	        (!boost_before_infection && sampling_time >= infection_times[x])){
          
          //      if(sampling_time >= infection_times[x]){
      	infection_time = infection_times[x];
          
       //   Rcpp::Rcout << std::endl << std::endl << "Infection time: " << infection_time << std::endl;
          
      	time = sampling_time - infection_time; // Time between sample and infection
      	inf_map_index = infection_strain_indices_tmp[x]; // Index of this infecting strain in antigenic map
      	
      	// Get titre for each strain from this infections
      	for(int k = 0; k < n_titres; ++k){
      	 // Rcpp::Rcout << "Strain index: " << k << std::endl;
           	// Add up contribution of all previous infections to titre that
      	    // would be observed at this infection time
      	  //  Rcpp::Rcout << "Going through previous infections..." << std::endl;
      	  	  monitored_titre = 0;
      	    for(int ii = x - 1; ii >= 0; --ii){
      	      
      	      // Get antigenic relationship between this infecting strain and the measured strain
      	      inf_map_index_tmp = measurement_strain_indices[tmp_titre_index + k]*number_strains + infection_strain_indices_tmp[ii];
      	      // inf_map_index_tmp = inf_map_index * number_strains + infection_strain_indices_tmp[ii];
      	      seniority = MAX(0, 1.0 - tau * ii);
      	      wane_amount = MAX(0, 1.0 - wane * (infection_time - infection_times[ii]));
      	      
      	     // Rcpp::Rcout << "This infection time: " << infection_times[ii] << std::endl;
      	     // Rcpp::Rcout << "Long term cr: " << antigenic_map_long[inf_map_index_tmp] << std::endl;
      	     // Rcpp::Rcout << "Titre at time of boost: " << monitored_titres(k,ii) << std::endl;
      	      
      	      long_boost = seniority * mu * antigenic_map_long[inf_map_index_tmp];
      	      short_boost = seniority * mu_short * antigenic_map_short[inf_map_index_tmp];
          	  if(monitored_titres(ii,k) >= boost_limit){
          	    long_boost *= titre_suppression;
          	    short_boost *= titre_suppression;
          	  } else {
          	    long_boost *= MAX(0,1.0 - gradient*monitored_titres(ii,k));
          	    short_boost *= MAX(0,1.0 - gradient*monitored_titres(ii,k));	    
          	  }
          	  long_boost = MAX(0, long_boost);
          	  short_boost = MAX(0, short_boost);
          	  boost = long_boost + short_boost * wane_amount;
          	  monitored_titre += boost;
      	    }
        	  monitored_titres(x,k) = monitored_titre;
      	}
      	wane_amount= MAX(0, 1.0 - (wane*time)); // Basic waning function
      	seniority = MAX(0, 1.0 - tau*n_inf); // Antigenic seniority
      	// Rcpp::Rcout << std::endl << std::endl << "Now boost from this infection..." << std::endl;
      	// Find contribution to each measured titre from this infection
      	for(int k = 0; k < n_titres; ++k){
      	  // Rcpp::Rcout << "Strain index: " << k << std::endl;
      	  index = measurement_strain_indices[tmp_titre_index + k]*number_strains + inf_map_index;
      
      	  long_boost = seniority * mu * antigenic_map_long[index];
      	  short_boost = seniority * mu_short * antigenic_map_short[index];
      	  
      	  // Rcpp::Rcout << "Long term cr: " << antigenic_map_long[index] << std::endl;
      	 // Rcpp::Rcout << "Titre at time of boost: " << monitored_titres(k,x) << std::endl;
      	//  Rcpp::Rcout << "Titre suppression: " << 1 - gradient * monitored_titres(k,x) << std::endl;
       // Titre dependent boosting - at ceiling
      	  if(monitored_titres(x,k) >= boost_limit){
      	    long_boost *= titre_suppression;
      	    short_boost *= titre_suppression;
      // Titre dependent boosting - below ceiling
      	  } else {
      	    long_boost = long_boost * (1 - gradient * monitored_titres(x,k)); 
      	    short_boost = short_boost * (1 - gradient * monitored_titres(x,k)); // Titre dependent boosting - below ceiling
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
 // Rcpp::Rcout << "End" << std::endl;
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


//' Titre dependent boosting fast
//' 
//' A fast implementation of the titre dependent boosting function, giving predicted titres for a number of samples for one individual. Note that this version attempts to minimise memory allocations.
//' @family boosting_functions
//' @seealso \code{\link{titre_data_fast}}
void titre_data_fast_individual_titredep_old(NumericVector &predicted_titres,
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
            long_boost = long_boost * MAX(0,1 - gradient * monitored_titres[x]); 
            short_boost = short_boost * MAX(0,1 - gradient * monitored_titres[x]); // Titre dependent boosting - below ceiling
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
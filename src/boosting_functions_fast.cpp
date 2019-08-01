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




//' Backboosting model
//' 
//' Back boosting model implementation
//' @family boosting_functions
//' @seealso \code{\link{titre_data_fast}}
void titre_model_backboost_cpp(NumericVector &predicted_titres,
			       const double &mu,
			       const double &mu_short,
			       const double &wane,
			       const double &tau,
			       const double &affinity_maturation,
			       const double &nu_long_recall,
			       const double &nu_short_recall,
			       const double &max_interference,
			       const double &interference_gradient,		      
			       const NumericVector &infection_times,
			       const IntegerVector &infection_strain_indices,
			       const IntegerVector &measurement_strain_indices,
			       const NumericVector &sample_times,
			       const int &index_in_samples,
			       const int &end_index_in_samples,
			       const int &start_index_in_data1,
			       const IntegerVector &nrows_per_blood_sample,
			       const int &number_strains,
			       const NumericVector &antigenic_map_short,
			       const NumericVector &antigenic_map_long,
			       const NumericVector &antigenic_distances,
			       bool boost_before_infection = false){

  // If mu is 0, then we're using nu as the main boost as well
  double use_mu = mu;
  double use_mu_short = mu_short;
  if(mu == 0) use_mu = nu_long_recall;
  if(mu_short == 0) use_mu_short = nu_long_recall;
  /*  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << "mu: " << mu << std::endl;
  Rcpp::Rcout << "mu_short: " << mu_short << std::endl;
  Rcpp::Rcout << "nu_long: " << nu_long_recall << std::endl;
  Rcpp::Rcout << "nu_short: " << nu_short_recall << std::endl;
  Rcpp::Rcout << "wane: " << wane << std::endl;
  Rcpp::Rcout << "max_interference: " << max_interference << std::endl;
  Rcpp::Rcout << "interference_gradient: " << interference_gradient << std::endl;
  Rcpp::Rcout << "tau: " << tau << std::endl;
  */
  int n_titres;
  int max_infections = infection_times.size();
  
  //NumericMatrix titres(n_titres, n_samps);
  
  // Index parameters
  int end_index_in_data;
  int tmp_titre_index;
  int start_index_in_data = start_index_in_data1;
  int inf_map_index; // Index of infection strain
  int index;
  
  // Time parameters
  double infection_time;
  double sampling_time;
  double time; // Time elapsed between sample and infection time
  double n_inf;
  
  // Waning and affinity maturation
  double wane_amount;
  double maturation_amount;
  
  // Parameters to calculate recall boosting amounts
  double overall_suppression;
  double distance;
  double boost_propn_tmp;
  double tmp_boost_long;
  double tmp_boost_short;
  double tmp_boost_propn;
  
  // Store boosting from novel and recall responses 
  NumericMatrix boost_propns(max_infections, max_infections);
  NumericMatrix long_term_boosts(max_infections, max_infections);
  NumericMatrix short_term_boosts(max_infections, max_infections);

  NumericVector seniority_terms(max_infections);
  double seniority;
  // Go through each infection and calculate HOMOLOGOUS boosting from novel and recall responses
  n_inf = 1.0;
  for(int i = 0; i < max_infections; ++i){
    overall_suppression = 0;
    seniority_terms[i] = MAX(0, 1.0 - tau*(n_inf - 1.0));
    n_inf++;
    // From this infection, calculate back-boost to each previous strain
    for(int j = 0; j < i; ++j){
      // Find distance between this infection and each previous strain
      distance = antigenic_distances[infection_strain_indices[i]*number_strains + infection_strain_indices[j]];
      boost_propn_tmp = MAX(0, max_interference - interference_gradient*distance);
      
      // Total boosting stolen
      overall_suppression += boost_propn_tmp;
      
      // If recall boost, use recall pars
      long_term_boosts(i, j) = nu_long_recall;
      short_term_boosts(i, j) = nu_short_recall; // Does recall have antigenically broad short term boost?
      boost_propns(i, j) = boost_propn_tmp;
    }
    // If there was some stolen recall boost, re-normalise
    if(overall_suppression > 0){
      boost_propns(i,_) = boost_propns(i,_)/overall_suppression;
      boost_propns(i,_) = boost_propns(i,_)*MIN(max_interference, overall_suppression);
    }
    //Rcpp::Rcout << "Back boosts: " << boost_propns(i,_) << std::endl;
    // Remaining boost is novel boost
    boost_propns(i,i) = 1.0 - MIN(max_interference, overall_suppression);
    long_term_boosts(i,i) = mu;
    short_term_boosts(i,i) = mu_short;
  }

  // Now that boosting has been calculated, work out observed titre for each observed strain at
  // each sample time this individual has
  // For each sample
  
  // l = recalled strain
  // j = current infecting strain
  // i = sampling time
  // k = strain being measured
  
  for(int i = index_in_samples; i <= end_index_in_samples; ++i){
    sampling_time = sample_times[i];
    //n_inf = 1.0;
    n_titres = nrows_per_blood_sample[i];
    end_index_in_data = start_index_in_data + n_titres;
    tmp_titre_index = start_index_in_data;
    //Rcpp::Rcout << "Sampling time: " << sampling_time << std::endl;
    // Contribution of all infections
    for(int j = 0; j < max_infections; ++j){
      infection_time = infection_times[j];
      if((boost_before_infection && sampling_time > infection_time) ||
	 (!boost_before_infection && sampling_time >= infection_time)){
	// ####################
	// ## HERE TO CHANGE INFECTION GIVING BOOST AT SAMPLE TIME RATHER THAN AFTER
	// ####################
        // Time since infection 
        time = sampling_time - infection_time;
        
        // How much waning and affinity maturation of this recall response
        // will have happened by this sample time
        wane_amount = MAX(0, 1.0 - wane*time);
        maturation_amount = MIN(1.0, affinity_maturation);
	seniority = seniority_terms[j];
        // Each infection boosts recall response
        for(int l = 0; l <= j; ++l){
          // Which strain does this recall infection correspond to?
          inf_map_index = infection_strain_indices[l];
          tmp_boost_long = long_term_boosts(j, l);
	  tmp_boost_short = short_term_boosts(j, l);
          tmp_boost_propn = boost_propns(j,l);
          // For each measured strain at each sample
          for(int k = 0; k < n_titres; ++k){
            index = measurement_strain_indices[tmp_titre_index + k]*number_strains + inf_map_index;
            predicted_titres[tmp_titre_index + k] +=  seniority*tmp_boost_propn* // Modified boost
              (tmp_boost_long*antigenic_map_long[index]*maturation_amount + // Long term
               tmp_boost_short*antigenic_map_short[index]*wane_amount); // Short term
          }
        }
      }
    }
    start_index_in_data = end_index_in_data;
  }
}


//' Age specific boosting fast
//' 
//' A fast implementation of the age mediated boosting function, giving predicted titres for a number of samples for one individual. Note that this version attempts to minimise memory allocations.
//' @family boosting_functions
//' @seealso \code{\link{titre_data_fast}}
void titre_data_fast_individual_age(NumericVector &predicted_titres,
				    const double &mu,
				    const double &mu_short,
				    const double &age_gradient,
				    const double &age_min_boost_propn,
				    const double &wane,
				    const double &tau,
				    const double &age,
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
  double age_at_inf;
  double age_reduction;
    
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
	time = sampling_time - infection_times[x]; // Time between sample and infection
	age_at_inf = infection_times[x] - age;
	age_reduction = MAX(age_min_boost_propn, 1.0 - age_gradient*age_at_inf);
	wane_amount= MAX(0, 1.0 - (wane*time)); // Basic waning function
	seniority = MAX(0, 1.0 - tau*(n_inf - 1.0)); // Antigenic seniority
	inf_map_index = infection_strain_indices_tmp[x]; // Index of this infecting strain in antigenic map

	// Find contribution to each measured titre from this infection
	for(int k = 0; k < n_titres; ++k){
	  index = measurement_strain_indices[tmp_titre_index + k]*number_strains + inf_map_index;
	  predicted_titres[tmp_titre_index + k] += seniority * age_reduction *
	    ((mu*antigenic_map_long[index]) + (mu_short*antigenic_map_short[index])*wane_amount);
	}
	++n_inf;
      }
    }
    start_index_in_data = end_index_in_data;
  }
}


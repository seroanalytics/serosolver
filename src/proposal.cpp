#include <RcppArmadilloExtensions/sample.h>
#include "boosting_functions.h"
#include "infection_model.h"
#include "likelihood_funcs.h"
#include "helpers.h"
// [[Rcpp::depends(RcppArmadillo)]]


//' Infection history gibbs proposal, fast
//'  Generates a new infection history matrix and corresponding individual likelihoods, using a gibbs sampler from the infection history prior. See \code{\link{infection_history_proposal_gibbs}}, as inputs are very similar.
//' @param theta NumericVector, the model parameters used to solve the model
//' @param infection_history_mat IntegerMatrix the matrix of 1s and 0s corresponding to individual infection histories
//' @param old_probs_1 NumericVector, the current likelihoods for each individual
//' @param sampled_indivs IntegerVector, indices of sampled individuals
//' @param n_years_samp int, for each individual, how many time periods to resample infections for?
//' @param age_mask IntegerVector, length of the number of individuals, with indices specifying first time period that an individual can be infected (indexed from 1, such that a value of 1 allows an individual to be infected in any time period)
//' @param strain_mask IntegerVector, length of the number of individuals, with indices specifying last time period that an individual can be infected (ie. last time a sample was taken)
//' @param n_alive IntegerVector, length of the number of time periods that an individual could be infected, giving the number of individual alive in each time period
//' @param swap_propn double, what proportion of proposals should be swap steps (ie. swap contents of two cells in infection_history rather than adding/removing infections)
//' @param swap_distance int, in a swap step, how many time steps either side of the chosen time period to swap with
//' @param alpha double, alpha parameter for beta prior on infection probability
//' @param beta double, beta parameter for beta prior on infection probability
//' @param circulation_times NumericVector, the times that each strain circulated
//' @param circulation_times_indices IntegerVector, indexing vector from 0:(number of strains-1)
//' @param sample_times NumericVector, the vector of real times that samples were taken
//' @param rows_per_indiv_in_samples IntegerVector, How many rows in titre data correspond to each individual, sample and repeat?
//' @param cum_nrows_per_individual_in_data IntegerVector, How many rows in the titre data correspond to each individual?
//' @param cum_nrows_per_individual_in_repeat_data IntegerVector, For the repeat data (ie. already calculated these titres), how many rows in the titre data correspond to each individual?
//' @param nrows_per_blood_sample IntegerVector, Split the sample times and runs for each individual
//' @param measurement_strain_indices IntegerVector, For each titre measurement, corresponding entry in antigenic map
//' @param antigenic_map_long NumericVector, the collapsed cross reactivity map for long term boosting, after multiplying by sigma1, see \code{\link{create_cross_reactivity_vector}}
//' @param antigenic_map_short NumericVector, the collapsed cross reactivity map for short term boosting, after multiplying by sigma2, see \code{\link{create_cross_reactivity_vector}}
//' @param data NumericVector, data for all individuals for the first instance of each calculated titre
//' @param repeat_data NumericVector, the repeat titre data for all individuals (ie. do not solve the same titres twice)
//' @param repeat_indices IntegerVector, which index in the main data vector does each entry in repeat_data correspond to ie. which calculated titre in predicted_titres should be used for each observation?
//' @param titre_shifts NumericVector, if length matches the length of \code{data}, adds these as measurement shifts to the predicted titres. If lengths do not match, is not used.
//' @param mus NumericVector, if length is greater than one, assumes that strain-specific boosting is used rather than a single boosting parameter
//' @param boosting_vec_indices IntegerVector, same length as circulation_times, giving the index in the vector \code{mus} that each entry should use as its boosting parameter.
//' @param solve_likelihood bool, if FALSE does not solve likelihood when calculating acceptance probability
//' @param temp double, temperature for parallel tempering MCMC
//' @return an R list, one entry with the matrix of 1s and 0s corresponding to the infection histories for all individuals and the other with the new corresponding likelihoods per individual
//' @export
//' @family infection_history_proposal
// [[Rcpp::export]]
List infection_history_proposal_gibbs_fast(const NumericVector &theta, // Model parameters
					   const IntegerMatrix &infection_history_mat,  // Current infection history
					   const NumericVector &old_probs_1,
					   const IntegerVector &sampled_indivs,
					   const IntegerVector &n_years_samp_vec,
					   const IntegerVector &age_mask, // Age mask
					   const IntegerVector &strain_mask, // Age mask
					   const IntegerVector &n_alive, // Number of individuals alive in each year
					   const double &swap_propn,
					   const int &swap_distance,
					   const double &alpha, // Alpha for prior
					   const double &beta, // Beta for prior
					   const NumericVector &circulation_times,
					   const IntegerVector &circulation_times_indices,
					   const NumericVector &sample_times,
					   const IntegerVector &rows_per_indiv_in_samples, // How many rows in unique sample times table correspond to each individual?
					   const IntegerVector &cum_nrows_per_individual_in_data, // How many rows in the titre data correspond to each individual?
					   const IntegerVector &cum_nrows_per_individual_in_repeat_data, // How many rows in the repeat titre data correspond to each individual?
					   const IntegerVector &nrows_per_blood_sample, // How many rows in the titre data table correspond to each unique individual + sample time + repeat?
					   const IntegerVector &measurement_strain_indices, // For each titre measurement, corresponding entry in antigenic map
					   const NumericVector &antigenic_map_long, 
					   const NumericVector &antigenic_map_short,
					   const NumericVector &data,
					   const NumericVector &repeat_data,
					   const IntegerVector &repeat_indices,
					   const NumericVector &titre_shifts,
					   const NumericVector &mus,
					   const IntegerVector &boosting_vec_indices,
					   const double temp=1,
					   bool solve_likelihood=true
					   ){
  // ########################################################################
  // Parameters to control indexing of data
  IntegerMatrix new_infection_history_mat(infection_history_mat); // Can this be avoided? Create a copy of the inf hist matrix
  int n_titres_total = data.size(); // How many titres are there in total?
  NumericVector predicted_titres(n_titres_total); // Vector to store predicted titres
  NumericVector old_probs = clone(old_probs_1); // Create a copy of the current old probs

  // These can be pre-computed
  int n_indivs = infection_history_mat.nrow();  // How many individuals are there in total?
  int number_strains = infection_history_mat.ncol(); // How many possible years are we interested in?
  int n_sampled = sampled_indivs.size(); // How many individuals are we actually investigating?

  // These indices allow us to step through the titre data vector
  // as if it were a matrix ie. number of rows for each individual
  // at a time
  int index_in_samples; // Index in sample times vector to point to
  int end_index_in_samples; // Index in sample times vector to end at
  int start_index_in_data; // Index in titre data to start at
  int end_index_in_data; // Index in titre data to end at

  int start_index_in_repeat_data; // Index in repeat titre data to start at
  int end_index_in_repeat_data; // Index in repeat titre data to end at

  int tmp_titre_index; // Index in titre data that we are currently at
  int inf_map_index; // Index in antigenic map of the infecting strain
  int index; // 
  int number_samples; // Tmp store number of blood samples for this individual
  int n_titres; // Tmp store of number of titres to predict for this sample
 
  double sampling_time; // Tmp store time that blood sample taken
  double time; // Tmp store time between sample and exposure

  IntegerVector new_infection_history(number_strains); // New proposed infection history
  IntegerVector infection_history(number_strains); // Old infection history
  LogicalVector indices;

  NumericVector infection_times; // Tmp store infection times for this infection history, combined with indices
  IntegerVector infection_strain_indices_tmp; // Tmp store which index in antigenic map these infection times relate to

  // Extract parameters first to avoid allocation costs
  double mu = theta["mu"];
  double mu_short = theta["mu_short"];
  double wane = theta["wane"];
  double tau = theta["tau"];
  double seniority;
  double n_inf;

  int wane_type = theta["wane_type"]; 
  bool alternative_wane_func = wane_type == 1;
  double kappa;
  double t_change;


  bool titre_dependent_boosting = theta["titre_dependent"] == 1;
  double gradient;
  double boost_limit;
  
  bool strain_dep_boost = false;

  if (alternative_wane_func){
    kappa = theta["kappa"];
    t_change = theta["t_change"];
  }  

  if (titre_dependent_boosting) {
    gradient = theta["gradient"];
    boost_limit = theta["boost_limit"];
  }
  
  if (mus.size() > 1) {
    strain_dep_boost = true;    
  }

  bool base_function = !(alternative_wane_func || titre_dependent_boosting || strain_dep_boost);
  // ########################################################################

  // ########################################################################
  IntegerVector samps; // Variable vector to sample from
  IntegerVector locs; // Vector of locations that were sampled
  // ########################################################################

  // ########################################################################
  int indiv; // Index of the individual under consideration
  int year; // Index of year being updated

  int n_samp_max; // Maximum number of years to sample for individual
  int n_samp_length; // Number of years that COULD be sampled for this individual
  int new_entry;
  int old_entry;
  int loc1, loc2, tmp; // Which indices are we looking at?
  int n_years_samp; // How many years to sample for this individual?
  int loc1_val_old, loc2_val_old;

  // ########################################################################

  // ########################################################################
  double m; // number of infections in a given year
  double n; // number alive in a particular year

  double m_1_new, m_1_old,m_2_new,m_2_old;
  double n_1, n_2;
  double prior_1_old, prior_2_old, prior_1_new,prior_2_new,prior_new,prior_old;

  double rand1; // Store a random number
  double ratio; // Store the gibbs ratio for 0 or 1 proposal

  double old_prob; // Likelihood of old number
  double new_prob; // Likelihood of new number
  double log_prob; // Likelihood ratio

  double lbeta_const = R::lbeta(alpha, beta);

  // For likelihood
  const double sd = theta["error"];
  const double den = sd*M_SQRT2;
  const double max_titre = theta["MAX_TITRE"];
  const double log_const = log(0.5);


  // Titre dependent boosting
  bool use_titre_shifts = false;
  if(titre_shifts.size() == n_titres_total) use_titre_shifts = true;
  // ########################################################################
  
  // ########################################################################
  // For each individual
  for(int i = 0; i < n_sampled; ++i){
    // Get index of individual under consideration and their current likelihood
    indiv = sampled_indivs[i]-1;
    old_prob = old_probs_1[indiv];  
    // Indexing for data upkeep
    index_in_samples = rows_per_indiv_in_samples[indiv];
    end_index_in_samples = rows_per_indiv_in_samples[indiv+1] - 1;
    number_samples = end_index_in_samples - index_in_samples;      

    n_years_samp = n_years_samp_vec[indiv]; // How many years are we intending to resample from for this individual?
    n_samp_length  = strain_mask[indiv] - age_mask[indiv]; // How many years maximum can we sample from?
    n_samp_max = std::min(n_years_samp, n_samp_length); // Use the smaller of these two numbers
    // Indexing for data upkeep
    start_index_in_data = cum_nrows_per_individual_in_data[indiv];
    end_index_in_data = cum_nrows_per_individual_in_data[indiv+1]-1;
    start_index_in_repeat_data = cum_nrows_per_individual_in_repeat_data[indiv];

    // Swap contents of a year for an individual
    if(R::runif(0,1) < swap_propn){
      new_infection_history = new_infection_history_mat(indiv,_);

      loc1 = floor(R::runif(0,n_samp_length+1)); // Choose a location from age_mask to strain_mask
      loc2 = loc1 + floor(R::runif(-swap_distance,swap_distance)); // Perturb +/- swap_distance
 
      // If we have gone too far left or right, reflect at the boundaries
      // Reflect at boundary
      /*if(loc2 < 0 || loc2 >= n_samp_length){
	loc2 = n_samp_length - (loc2 - floor(loc2/n_samp_length)*n_samp_length) - 1;
	}*/
      
      
      while(loc2 < 0) loc2 += n_samp_length;
      if(loc2 >= n_samp_length) loc2 -= floor(loc2/n_samp_length)*n_samp_length;
      
      // Get onto right scale (starting at age mask)
      loc1 += age_mask[indiv] - 1;
      loc2 += age_mask[indiv] - 1;

      loc1_val_old = new_infection_history(loc1);
      loc2_val_old = new_infection_history(loc2);

      // Only proceed if we've actually made a change
      if(loc1_val_old != loc2_val_old){
	// Swap contents
	new_infection_history(loc1) = new_infection_history(loc2);
	new_infection_history(loc2) = loc1_val_old;

	// Prior for previous state
	m_1_old = sum(new_infection_history_mat(_,loc1));      
	m_2_old = sum(new_infection_history_mat(_,loc2));
	    
	// Number alive is number alive overall
	n_1 = n_alive(loc1);
	n_2 = n_alive(loc2);

	// Pre-compute these? 
	prior_1_old = R::lbeta(m_1_old + alpha, n_1 - m_1_old + beta)-lbeta_const;
	prior_2_old = R::lbeta(m_2_old + alpha, n_2 - m_2_old + beta)-lbeta_const;
	prior_old = prior_1_old + prior_2_old;

	// Prior for new state
	m_1_new = m_1_old - loc1_val_old + loc2_val_old;
	m_2_new = m_2_old - loc2_val_old + loc1_val_old;

	prior_1_new = R::lbeta(m_1_new + alpha, n_1 - m_1_new + beta)-lbeta_const;
	prior_2_new = R::lbeta(m_2_new + alpha, n_2 - m_2_new + beta)-lbeta_const;
	prior_new = prior_1_new + prior_2_new;	 

	////////////////////////
	if(solve_likelihood){
	  // Calculate likelihood!
	  indices = new_infection_history > 0;
	  infection_times = circulation_times[indices];
	  infection_strain_indices_tmp = circulation_times_indices[indices];
	  
	  if (base_function) {
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
						antigenic_map_long);	
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
							antigenic_map_long);
	  } else {
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
					     antigenic_map_long);
	  }

	  if(use_titre_shifts){
	    add_measurement_shifts(predicted_titres, titre_shifts, 
				   start_index_in_data, end_index_in_data);
	  }
	  
	  // Now have all predicted titres for this individual calculated
	  // Need to calculate likelihood of these titres... 
	  new_prob = 0;

	  // NEED SOMETHING TO DEAL WITH REPEATS
	  // Go from first row in the data for this individual to up to the next one, accumlating
	  // likelihood for this individual
	  // For unique data
	  proposal_likelihood_func(new_prob, predicted_titres, indiv, data, repeat_data, repeat_indices,
				   cum_nrows_per_individual_in_data, cum_nrows_per_individual_in_repeat_data,
				   log_const, den, max_titre);
	} else {
	  old_prob = new_prob = old_probs[indiv];
	}

	log_prob = std::min<double>(0.0, (new_prob+prior_new) - (old_prob+prior_old));
	rand1 = R::runif(0,1);
	if(log(rand1) < log_prob/temp){
	  // Update the entry in the new matrix Z
	  old_prob = new_prob;
	  old_probs[indiv] = new_prob;
	  // Carry out the swap
	  tmp = new_infection_history_mat(indiv,loc1);
	  new_infection_history_mat(indiv,loc1) = new_infection_history_mat(indiv,loc2);
	  new_infection_history_mat(indiv,loc2) = tmp;
	}
      }
    } else {
      // Sample n_samp_real years from 0:length. Ths will be used to pull years from
      // sample_years
      // Note sampling indices in the individual's infection history, not the matrix Z
      // Take n_samp_max random samples 
      samps = seq(0, n_samp_length);    // Create vector from 0:length of alive years
      locs = RcppArmadillo::sample(samps, n_samp_max, FALSE, NumericVector::create());

      for(int j = 0; j < n_samp_max; ++j){
	new_infection_history = new_infection_history_mat(indiv,_);
	year = locs[j] + age_mask[indiv] - 1;

	old_entry = new_infection_history(year);

	// Get number of individuals that were alive and/or infected in that year,
	// less the current individual
	// Number of infections in this year, less infection status of this individual in this year
	m = sum(new_infection_history_mat(_,year)) - old_entry;
	n = n_alive(year) - 1;
	       
	// Work out proposal ratio - prior from alpha, beta and number of other infections
	ratio = (m + alpha)/(n + alpha + beta);

	// Propose 1 or 0 based on this ratio
	rand1 = R::runif(0,1);

	if(rand1 < ratio){
	  new_entry = 1;
	  new_infection_history(year) = 1;
	} else {
	  new_entry = 0;
	  new_infection_history(year) = 0;
	}
	// If proposing a change, need to check likelihood ratio
	if(new_entry != old_entry){
	  if(solve_likelihood){
	    // Calculate likelihood!
	    indices = new_infection_history > 0;
	    infection_times = circulation_times[indices];
	    if(infection_times.size() > 0){
	      infection_strain_indices_tmp = circulation_times_indices[indices];
	      if (base_function) {
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
						    antigenic_map_long);
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
							    antigenic_map_long);	
	      } else {
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
						 antigenic_map_long);
	      }

	    }
	    
	    if(use_titre_shifts){
	      add_measurement_shifts(predicted_titres, titre_shifts, 
				     start_index_in_data, end_index_in_data);
	    }
	    
	    // Now have all predicted titres for this individual calculated
	    // Need to calculate likelihood of these titres... 
	    new_prob = 0;

	    // NEED SOMETHING TO DEAL WITH REPEATS
	    // Go from first row in the data for this individual to up to the next one, accumlating
	    // likelihood for this individual
	    // For unique data
	    proposal_likelihood_func(new_prob, predicted_titres, indiv, data, repeat_data, repeat_indices,
				     cum_nrows_per_individual_in_data, cum_nrows_per_individual_in_repeat_data,
				     log_const, den, max_titre);
	    // =====================
	  } else {
	    old_prob = new_prob = old_probs[indiv];
	  }
	  // Don't need to take into account prior prob, as sampling from this
	  log_prob = std::min<double>(0.0, new_prob - old_prob);
	  rand1 = R::runif(0,1);
	  if(log(rand1) < log_prob/temp){
	    // Update the entry in the new matrix Z
	    old_prob = new_prob;
	    old_probs[indiv] = new_prob;
	    new_infection_history_mat(indiv, year) = new_entry;
	  }
	}
      }
    }
  }
  List ret;
  ret["old_probs"] = old_probs;
  ret["new_infection_history"] = new_infection_history_mat;
  return(ret);
}




//' Gibbs sampling of infection histories
//'
//' Proposes a new matrix of infection histories by sampling from the prior on an individual's infection presence/absence in a particular time period, conditional on all other individuals in that time period. This allows us to integrate out the infection probability term for each time period. Should look at \code{\link{create_posterior_func}} for more details about the input parameters.
//' @param pars NumericVector, the model parameters used to solve the model, extract likelihood function parameters and alpha/beta
//' @param infection_history_mat IntegerMatrix the matrix of 1s and 0s corresponding to individual infection histories
//' @param indiv_samp_propn double, what proportion of individuals to resample in this proposal step
//' @param n_years_samp int, for each individual, how many time periods to resample infections for?
//' @param age_mask IntegerVector, length of the number of individuals, with indices specifying first time period that an individual can be infected (indexed from 1, such that a value of 1 allows an individual to be infected in any time period)
//' @param strain_mask IntegerVector, length of the number of individuals, with indices specifying last time period that an individual can be infected (ie. last time a sample was taken)
//' @param n_alive IntegerVector, length of the number of time periods that an individual could be infected, giving the number of individual alive in each time period
//' @param swap_propn double, what proportion of proposals should be swap steps (ie. swap contents of two cells in infection_history rather than adding/removing infections)
//' @param swap_distance int, in a swap step, how many time steps either side of the chosen time period to swap with
//' @param alpha double, alpha parameter for beta prior on infection probability
//' @param beta double, beta parameter for beta prior on infection probability
//' @param circulation_times NumericVector, the times that each strain circulated
//' @param circulation_times_indices IntegerVector, indexing vector from 0:(number of strains-1)
//' @param sample_times NumericVector, the vector of real times that samples were taken
//' @param rows_per_indiv_in_samples IntegerVector, How many rows in titre data correspond to each individual, sample and repeat?
//' @param cum_nrows_per_individual_in_data IntegerVector, How many rows in the titre data correspond to each individual?
//' @param nrows_per_blood_sample IntegerVector, Split the sample times and runs for each individual
//' @param measurement_strain_indices IntegerVector, For each titre measurement, corresponding entry in antigenic map
//' @param antigenic_map_long NumericVector, the collapsed cross reactivity map for long term boosting, after multiplying by sigma1
//' @param antigenic_map_short NumericVector, the collapsed cross reactivity map for short term boosting, after multiplying by sigma2
//' @param data NumericVector, the titre data for all individuals
//' @param to_add Nullable<NumericVector>, optional vector of measurement shifts to apply to all titres
//' @param additional_arguments Nullable<List>, optional list to pass down to titre solving function
//' @param DOBs NumericVector, vector of ages for each individual
//' @param solve_likelihood bool, if FALSE does not solve likelihood when calculating acceptance probability
//' @param total_alive int, if a positive number, uses this rather than the vector of alive individuals for the infection history prior
//' @param temp double, temperature for parallel tempering MCMC
//' @return a matrix of 1s and 0s corresponding to the infection histories for all individuals
//' @family infection_history_proposal
// [[Rcpp::export]]
IntegerMatrix infection_history_proposal_gibbs(const NumericVector& pars,
					       const IntegerMatrix& infection_history_mat,
					       double indiv_samp_propn,
					       const IntegerVector& n_years_samp_vec,
 					       const IntegerVector& age_mask,
 					       const IntegerVector& strain_mask,
					       const IntegerVector& n_alive,
					       double swap_propn,
					       int swap_distance,
					       double alpha,
					       double beta,
					       const NumericVector &circulation_times,
					       const IntegerVector &circulation_times_indices,
					       const NumericVector &sample_times,
					       const IntegerVector &rows_per_indiv_in_samples,
					       const IntegerVector &cum_nrows_per_individual_in_data,
					       const IntegerVector &nrows_per_blood_sample,
					       const IntegerVector &measurement_strain_indices,
					       const NumericVector &antigenic_map_long, 
					       const NumericVector &antigenic_map_short,
					       const NumericVector &data,
					       const Nullable<NumericVector> &to_add,
					       const Nullable<List> &additional_arguments,
					       const NumericVector &DOBs,
					       bool solve_likelihood=true,
					       int total_alive=-1,
					       double temp=1
					       ){
  // ########################################################################
  // Parameters to control indexing of data
  IntegerMatrix new_infection_history_mat(infection_history_mat);
  int n_indivs = new_infection_history_mat.nrow();
  int n_strains = new_infection_history_mat.ncol();
  // These indices allow us to step through the titre data vector
  // as if it were a matrix ie. number of rows for each individual
  // at a time
  int start_index_samples;
  int end_index_samples;
  int start_index_data;
  int end_index_data;
  // ########################################################################

  // ########################################################################
  IntegerVector infs_in_year;
  IntegerVector samps;
  IntegerVector sample_years;
  IntegerVector locs;
  IntegerVector proposed_indiv_hist;
  IntegerVector indiv_hist;
  IntegerVector swap_size;
  // ########################################################################

  // ########################################################################
  int indiv; // Index of the individual under consideration
  int year; // Index of year being updated
  double DOB; // The DOB of the individual for solving the model, if needed
  int n_samp_max; // Maximum number of years to sample
  int new_entry;
  int old_entry;
  int n_infected=0;
  int loc1, loc2, tmp;
  int n_years_samp;
  // ########################################################################

  // ########################################################################
  double m; // number of infections in a given year
  double n; // number alive in a particular year

  //double alpha;
  //double beta;

  double m_1_new, m_1_old,m_2_new,m_2_old;
  double n_1_new, n_1_old, n_2_new, n_2_old;
  double prior_1_old, prior_2_old, prior_1_new,prior_2_new,prior_new,prior_old;

  double rand1; // Store a random number
  double ratio; // Store the gibbs ratio for 0 or 1 proposal

  double old_prob; // Likelihood of old number
  double new_prob; // Likelihood of new number
  double log_prob; // Likelihood ratio
  // ########################################################################
  
  bool prior_on_total = total_alive > 0;
  // ########################################################################
  NumericVector to_add_tmp(1);
  //NumericVector to_add_tmp;

  if(prior_on_total){
    n_infected = sum(infection_history_mat);
  }
  // ########################################################################
  // For each individual
  for(int i = 1; i <= n_indivs; ++i){
    // Indexing for data upkeep
    start_index_samples = rows_per_indiv_in_samples[i-1];
    end_index_samples = rows_per_indiv_in_samples[i] - 1;

    start_index_data = cum_nrows_per_individual_in_data[i-1];
    end_index_data = cum_nrows_per_individual_in_data[i] - 1;
    if (to_add.isNotNull()){
      to_add_tmp = subset_nullable_vector(to_add, start_index_data, end_index_data);
    }
    
    // Choose whether to sample this individual or not
    if(R::runif(0,1) < indiv_samp_propn){
      // Index of this individual
      indiv = i-1;
      DOB = DOBs[indiv];
      n_years_samp = n_years_samp_vec[indiv];
      // Make vector of year indices to sample from
      // These are the indices in the matrix Z
      sample_years = seq(age_mask[indiv]-1,strain_mask[indiv]-1);
      // Sample the minimum of either the number of years alive, or 
      // the number of years that are asked to be changed
      n_samp_max = sample_years.size();
      samps = seq(0, n_samp_max-1);    // Create vector from 0:length of alive years
      n_samp_max = std::min(n_years_samp, n_samp_max);
      // Swap contents of a year for an individual
      if(R::runif(0,1) < swap_propn){
	proposed_indiv_hist = new_infection_history_mat(indiv,_);
	indiv_hist = new_infection_history_mat(indiv,_);
	  
	loc1 = samps(floor(R::runif(0,samps.size())));
	swap_size = seq(-swap_distance,swap_distance);
	loc2 = loc1 + swap_size(floor(R::runif(0,swap_size.size())));

	while(loc2 < 0) loc2 = loc2 + samps.size();
	if(loc2 >= samps.size()) loc2 = loc2 - floor(loc2/samps.size())*samps.size();
	loc1 = sample_years(loc1);
	loc2 = sample_years(loc2);

	if(new_infection_history_mat(indiv,loc1) != new_infection_history_mat(indiv, loc2)){
	  tmp = proposed_indiv_hist(loc1);
	  proposed_indiv_hist(loc1) = proposed_indiv_hist(loc2);
	  proposed_indiv_hist(loc2) = tmp;

	  if(!prior_on_total){
	    // Prior for previous state
	    m_1_old = sum(new_infection_history_mat(_,loc1));      
	    m_2_old = sum(new_infection_history_mat(_,loc2));
	    
	    // Number alive is number alive overall
	    n_1_old = n_alive(loc1);
	    n_2_old = n_alive(loc2);

	    prior_1_old = R::lbeta(m_1_old + alpha, n_1_old - m_1_old + beta)-R::lbeta(alpha,beta);
	    prior_2_old = R::lbeta(m_2_old + alpha, n_2_old - m_2_old + beta)-R::lbeta(alpha,beta);
	    prior_old = prior_1_old + prior_2_old;

	    // Prior for new state
	    m_1_new = sum(new_infection_history_mat(_,loc1)) - new_infection_history_mat(indiv,loc1) + proposed_indiv_hist(loc1);      
	    m_2_new = sum(new_infection_history_mat(_,loc2)) - new_infection_history_mat(indiv,loc2) + proposed_indiv_hist(loc2);

	    // Number alive is number alive overall
	    n_1_new = n_alive(loc1);
	    n_2_new = n_alive(loc2);

	    prior_1_new = R::lbeta(m_1_new + alpha, n_1_new - m_1_new + beta)-R::lbeta(alpha,beta);
	    prior_2_new = R::lbeta(m_2_new + alpha, n_2_new - m_2_new + beta)-R::lbeta(alpha,beta);
	    prior_new = prior_1_new + prior_2_new;
	  } else {
	    prior_old = prior_new = 0;
	  }

	  if(solve_likelihood){
	    old_prob = likelihood_data_individual(pars, 
						  indiv_hist, 
						  circulation_times, 
						  circulation_times_indices,
						  sample_times[Range(start_index_samples, end_index_samples)], 
						  nrows_per_blood_sample[Range(start_index_samples,end_index_samples)], 
						  measurement_strain_indices[Range(start_index_data,end_index_data)],
						  antigenic_map_long, 
						  antigenic_map_short,  
						  n_strains,
						  data[Range(start_index_data,end_index_data)],
						  to_add_tmp,
						  DOB,
						  additional_arguments);

	    new_prob = likelihood_data_individual(pars, proposed_indiv_hist,
						  circulation_times, circulation_times_indices,
						  sample_times[Range(start_index_samples, end_index_samples)], 
						  nrows_per_blood_sample[Range(start_index_samples,end_index_samples)], 
						  measurement_strain_indices[Range(start_index_data,end_index_data)],
						  antigenic_map_long, 
						  antigenic_map_short,  
						  n_strains,
						  data[Range(start_index_data,end_index_data)],
						  to_add_tmp,
						  DOB,
						  additional_arguments);
	  } else {
	    old_prob = new_prob = 0;
	  }

	  log_prob = std::min<double>(0.0, (new_prob+prior_new)/temp - (old_prob+prior_old));
	  rand1 = R::runif(0,1);
	  if(log(rand1) < log_prob){
	    // Update the entry in the new matrix Z
	    tmp = new_infection_history_mat(indiv,loc1);
	    new_infection_history_mat(indiv,loc1) = new_infection_history_mat(indiv,loc2);
	    new_infection_history_mat(indiv,loc2) = tmp;
	  }
	}
      } else {
	// Sample n_samp_real years from 0:length. Ths will be used to pull years from
	// sample_years
	// Note sampling indices in the individual's infection history, not the matrix Z
	locs = RcppArmadillo::sample(samps, n_samp_max, FALSE, NumericVector::create());
	// ===========
	// This bit of code can be reactivated to propose multiple infection history changes at once
	// ===========
	// For each sampled year
	/*int k = sum(new_infection_history_mat);
	  int nn = total_alive;
	  for(int ii = 0; ii < locs.size(); ++ii){
	  k -= new_infection_history_mat(indiv, sample_years(locs[ii]));
	  }
	  nn -= locs.size();
	  proposed_indiv_hist = new_infection_history_mat(indiv,_);
	  indiv_hist = new_infection_history_mat(indiv,_);
	*/
	for(int j = 0; j < n_samp_max; ++j){
	  proposed_indiv_hist = new_infection_history_mat(indiv,_);
	  indiv_hist = new_infection_history_mat(indiv,_);
	
	  // If using multiple infection state proposals
	  //  ratio = (k + alpha)/(nn + alpha + beta);

	  // Which year under consideration
	  // Get index for the matrix Z (not index for individual's inf hist)
	  // Note that locs[j] is the index in the individual's inf hist
	  year = sample_years(locs[j]);

	  // If we want a different alpha/beta per year
	  // alpha = alphas[year];
	  // beta = betas[year];

	  if(prior_on_total){
	    old_entry = new_infection_history_mat(indiv, year);
	    // NUmber infected overall
	    m = n_infected - old_entry;
	    // Number alive is number alive overall
	    n = total_alive - 1;
	  } else {
	    // Get number of individuals that were alive and/or infected in that year,
	    // less the current individual
	    // Number of infections in this year, less infection status of this individual in this year
	    m = sum(new_infection_history_mat(_,year)) - new_infection_history_mat(indiv,year);     
	    n = n_alive(year) - 1;
	  }
	  // Work out proposal ratio - prior from alpha, beta and number of other infections
	  ratio = (m + alpha)/(n + alpha + beta);
	  // Propose 1 or 0 based on this ratio
	  rand1 = R::runif(0,1);
	  if(rand1 < ratio){
	    new_entry = 1;
	    //k++;
	    proposed_indiv_hist(year) = 1;
	  } else {
	    new_entry = 0;
	    proposed_indiv_hist(year) = 0;
	  }
	  //nn++;
	
	  // If proposing a change, need to check likelihood ratio
	  if(new_entry != new_infection_history_mat(indiv,year)){
	    if(solve_likelihood){
	      old_prob = likelihood_data_individual(pars, indiv_hist, circulation_times, circulation_times_indices,
						    sample_times[Range(start_index_samples, end_index_samples)], 
						    nrows_per_blood_sample[Range(start_index_samples,end_index_samples)], 
						    measurement_strain_indices[Range(start_index_data,end_index_data)],
						    antigenic_map_long, 
						    antigenic_map_short,  
						    n_strains,
						    data[Range(start_index_data,end_index_data)],
						    to_add_tmp,
						    DOB,
						    additional_arguments);
	      new_prob = likelihood_data_individual(pars, proposed_indiv_hist,
						    circulation_times, circulation_times_indices,
						    sample_times[Range(start_index_samples, end_index_samples)], 
						    nrows_per_blood_sample[Range(start_index_samples,end_index_samples)], 
						    measurement_strain_indices[Range(start_index_data,end_index_data)],
						    antigenic_map_long, 
						    antigenic_map_short,  
						    n_strains,
						    data[Range(start_index_data,end_index_data)],
						    to_add_tmp,
						    DOB,
						    additional_arguments);
	    } else {
	      old_prob = new_prob = 0;
	    }
	    log_prob = std::min<double>(0.0, new_prob - old_prob);
	    rand1 = R::runif(0,1);
	    if(log(rand1) < log_prob/temp){
	      if(prior_on_total){
		n_infected = m + new_entry;
	      }
	      // Update the entry in the new matrix Z
	      new_infection_history_mat(indiv, year) = new_entry;
	    }
	  }
	}
      }
    }
  }
  return(new_infection_history_mat);
}

//' Fast infection history proposal function
//' 
//' Proposes a new matrix of infection histories using a beta binomial proposal distribution. This particular implementation allows for n_infs epoch times to be changed with each function call. Furthermore, the size of the swap step is specified for each individual by move_sizes.
//' @param infection_history_mat and RcppArmadillo matrix of infection histories, where rows represent individuals and columns represent potential infection times. The contents should be a set of 1s (presence of infection) and 0s (absence of infection)
//' @param sampled_indivs IntegerVector, indices of which individuals to resample. Note that this is indexed from 1 (ie. as if passing straight from R)
//' @param age_mask IntegerVector, for each individual gives the first column in the infection history matrix that an individual could have been exposed to indexed from 1. ie. if alive for the whole period, entry would be 1. If alive for the 11th epoch, entry would be 11.
//' @param move_sizes IntegerVector, how far can a swap step sample from specified for each individual
//' @param n_infs IntegerVector, how many infections to add/remove/swap with each proposal step for each individual
//' @param alpha double, alpha parameter of the beta binomial
//' @param beta double, beta parameter of the beta binomial
//' @param rand_ns NumericVector, a vector of random numbers for each sampled individual. The idea is to pre-specify whether an individual experiences an add/remove step or a swap step to avoid random number sampling in C++
//' @return a matrix of 1s and 0s corresponding to the infection histories for all individuals
//' @export
//' @family infection_history_proposal
// [[Rcpp::export]]
arma::mat inf_hist_prop_cpp(arma::mat infection_history_mat, 
			    const IntegerVector& sampled_indivs, 
			    const IntegerVector& age_mask,
			    const IntegerVector& strain_mask,
			    const IntegerVector& move_sizes, 
			    const IntegerVector& n_infs,
			    double alpha, 
			    double beta, 
			    const NumericVector& rand_ns) {
  // Copy input matrix
  arma::mat new_infection_history_mat = infection_history_mat;
  IntegerVector locs; // Locations to be updated
  arma::uvec locs1;
  arma::mat x;
  arma::mat y;
  IntegerVector samps;
  int max_i_indiv;
  int indiv;
  int k;
  int n_inf;
  int n;
  int move_max;
  int move;
  int id1;
  int id2;
  int tmp;
  int n_samp_max;

  double rand1;
  double ratio;
  
  // For each sampled individual
  for(int i = 0; i < sampled_indivs.size(); ++i){

    // Isolate that individual's infection histories
    indiv = sampled_indivs[i]-1;
    n_inf = n_infs[indiv];
    x = new_infection_history_mat.submat(indiv, age_mask[indiv]-1, indiv, strain_mask[indiv]-1);
    samps = seq_len(x.n_cols);
    n_samp_max = samps.size();
    n_samp_max = std::min(n_inf, n_samp_max);
    
    // With 50% probability, add/remove infections or swap infections
    if(rand_ns[i] < 1.0/2.0){
      // Sample N random locations
      locs = RcppArmadillo::sample(samps, n_samp_max, FALSE, NumericVector::create());
      locs1 = as<arma::uvec>(locs)-1;
      y = x.elem(locs1);
      // Count the number of 1s and 0s
      k = accu(x) - accu(y);
      n = x.size() - n_samp_max;
      
      // For each sampled location, choose to turn into a 1 or 0 depending
      // on the beta binomial distribution.
      for(int j = 0; j < n_samp_max; ++j){
        ratio = (alpha + k)/(alpha + beta + n);
        rand1 = R::runif(0,1);
	// With probability 'ratio', add a 1. ie. if many 1s already, less likely
	// to add more 1s depending on alpha and beta
        if(rand1 < ratio){
          x(locs1(j)) = 1;
          k++;
        } else {
          x(locs1(j)) = 0;
        }
        n++;
      }
    } else {
      // Otherwise, swap the contents of N random locations
      max_i_indiv = x.size();
      for(int j = 0; j < n_samp_max; j++){
        id1 = floor(R::runif(0,1)*x.size());
        move_max = move_sizes[indiv];
        move = floor(R::runif(0,1)*2*move_max) - move_max;
        id2 = id1 + move;
	while(id2 < 0) id2 += max_i_indiv;
	if(id2 >= max_i_indiv) id2 = (id2 % max_i_indiv);
	tmp = x[id1];
        x[id1] = x[id2];
        x[id2] = tmp;
      }
    }
    new_infection_history_mat.submat(indiv, age_mask[indiv]-1, indiv,  strain_mask[indiv]-1) = x;
  }
  return(new_infection_history_mat);
}


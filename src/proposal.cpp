#include <RcppArmadilloExtensions/sample.h>
#include "boosting_functions_fast.h"
#include "likelihood_funcs.h"
#include "helpers.h"
// [[Rcpp::depends(RcppArmadillo)]]

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
arma::mat inf_hist_prop_prior_v3(arma::mat infection_history_mat, 
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




//' Infection history gibbs proposal, fast
//'  Generates a new infection history matrix and corresponding individual likelihoods, using a gibbs sampler from the infection history prior. See \code{\link{infection_history_proposal_gibbs}}, as inputs are very similar.
//' @param theta NumericVector, the model parameters used to solve the model
//' @param infection_history_mat IntegerMatrix the matrix of 1s and 0s corresponding to individual infection histories
//' @param old_probs_1 NumericVector, the current likelihoods for each individual
//' @param sampled_indivs IntegerVector, indices of sampled individuals
//' @param n_years_samp int, for each individual, how many time periods to resample infections for?
//' @param age_mask IntegerVector, length of the number of individuals, with indices specifying first time period that an individual can be infected (indexed from 1, such that a value of 1 allows an individual to be infected in any time period)
//' @param strain_mask IntegerVector, length of the number of individuals, with indices specifying last time period that an individual can be infected (ie. last time a sample was taken)
//' @param n_alive IntegerMatrix, number of columns is the number of time periods that an individual could be infected, giving the number of individual alive in each time period. Number of rows is the number of distinct groups.
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
//' @param group_id_vec IntegerVector, vector with 1 entry per individual, giving the group ID of that individual
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
List inf_hist_prop_prior_v2_and_v4(const NumericVector &theta, // Model parameters
				   const IntegerMatrix &infection_history_mat,  // Current infection history
				   const NumericVector &old_probs_1,
				   const IntegerVector &sampled_indivs,
				   const IntegerVector &n_years_samp_vec,
				   const IntegerVector &age_mask, // Age mask
				   const IntegerVector &strain_mask, // Age mask
				   const IntegerMatrix &n_alive, // No. of individuals alive each year/group
				   IntegerMatrix &n_infections, // No. of infections in each year/group
				   IntegerVector &n_infected_group,
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
				   const IntegerVector &group_id_vec, // Which group does each individual belong to?
				   const IntegerVector &measurement_strain_indices, // For each titre measurement, corresponding entry in antigenic map
				   const NumericVector &antigenic_map_long, 
				   const NumericVector &antigenic_map_short,
				   const NumericVector &data,
				   const NumericVector &repeat_data,
				   const IntegerVector &repeat_indices,
				   const NumericVector &titre_shifts,
				   const NumericVector &mus,
				   const IntegerVector &boosting_vec_indices,
				   const IntegerVector &total_alive,
				   const double temp=1,
				   bool solve_likelihood=true				   
				   ){
  // ########################################################################
  // Parameters to control indexing of data
  IntegerMatrix new_infection_history_mat(infection_history_mat); // Can this be avoided? Create a copy of the inf hist matrix
  int n_titres_total = data.size(); // How many titres are there in total?
  NumericVector predicted_titres(n_titres_total); // Vector to store predicted titres
  NumericVector old_probs = clone(old_probs_1); // Create a copy of the current old probs

  // Variables related to solving likelihood and model as little as possible
  bool swap_step_option = true;
  bool lik_changed = true;
  
  // These quantities can be pre-computed
  int n_indivs = infection_history_mat.nrow();  // How many individuals are there in total?
  int number_strains = infection_history_mat.ncol(); // How many possible years are we interested in?
  int n_sampled = sampled_indivs.size(); // How many individuals are we actually investigating?
  int n_infected = 0; // How many individuals were infected?
  
  // Using prior version 2 or 4?
  bool prior_on_total = total_alive(0) > 0;
    
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

  int group_id; // Vector of group IDs for each individual
  
  double sampling_time; // Tmp store time that blood sample taken
  double time; // Tmp store time between sample and exposure

  IntegerVector new_infection_history(number_strains); // New proposed infection history
  IntegerVector infection_history(number_strains); // Old infection history
  LogicalVector indices;

  NumericVector infection_times; // Tmp store infection times for this infection history, combined with indices
  IntegerVector infection_strain_indices_tmp; // Tmp store which index in antigenic map these infection times relate to

   
  // ########################################################################
  // Parameters related to infection history sampling
  int indiv; // Index of the individual under consideration
  IntegerVector samps; // Variable vector to sample from
  IntegerVector locs; // Vector of locations that were sampled
  int year; // Index of year being updated
  int n_samp_max; // Maximum number of years to sample for individual
  int n_samp_length; // Number of years that COULD be sampled for this individual
  int new_entries;
  int old_entries;
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
  
  // ====================================================== //
  // =============== SETUP MODEL PARAMETERS =============== //
  // ====================================================== //
  // 1. Extract general parameters that apply to all models
  // Pull out model parameters so only need to allocate once
  double mu = theta["mu"];
  double mu_short = theta["mu_short"];
  double wane = theta["wane"];
  double tau = theta["tau"];
  double seniority;
  double n_inf;

  // 2. Extract model parameters that are for specific mechanisms
  //    set a boolean flag to choose between model versions
  
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
  
  // 3. If not using one of the specific mechanism functions, set the base_function flag to TRUE
  bool base_function = !(alternative_wane_func || titre_dependent_boosting || strain_dep_boost);
  
  // 4. Extra titre shifts
  bool use_titre_shifts = false;
  if(titre_shifts.size() == n_titres_total) use_titre_shifts = true;
  // ########################################################################
  
  // ########################################################################
  // For each individual
  for(int i = 0; i < n_sampled; ++i){
    // Which proposal step to take and do we need to calculate the likelihood    
    swap_step_option = R::runif(0,1) < swap_propn;
    
    // Get index, group and current likelihood of individual under consideration
    indiv = sampled_indivs[i]-1;
    group_id = group_id_vec[indiv];
    old_prob = old_probs_1[indiv];

    // Indexing for data upkeep
    index_in_samples = rows_per_indiv_in_samples[indiv];
    end_index_in_samples = rows_per_indiv_in_samples[indiv+1] - 1;
    number_samples = end_index_in_samples - index_in_samples;      

    start_index_in_data = cum_nrows_per_individual_in_data[indiv];
    end_index_in_data = cum_nrows_per_individual_in_data[indiv+1]-1;
    start_index_in_repeat_data = cum_nrows_per_individual_in_repeat_data[indiv];
    
    // Time sampling control
    n_years_samp = n_years_samp_vec[indiv]; // How many times are we intending to resample for this individual?
    n_samp_length  = strain_mask[indiv] - age_mask[indiv]; // How many times maximum can we sample from?

    // If swap step, only doing one proposal for this individual
    if(swap_step_option){
      n_samp_max = 1;
    } else {
      n_samp_max = std::min(n_years_samp, n_samp_length); // Use the smaller of these two numbers
    }
    // Sample n_samp_real years from 0:length. Ths will be used to pull years from
    // sample_years
    // Note sampling indices in the individual's infection history, not the matrix Z
    // Take n_samp_max random samples 
    samps = seq(0, n_samp_length);    // Create vector from 0:length of alive years
    locs = RcppArmadillo::sample(samps, n_samp_max, FALSE, NumericVector::create());
    
    // For each selected infection history entry
    for(int j = 0; j < n_samp_max; ++j){
      // Assume that proposal hasn't changed likelihood until shown otherwise
      lik_changed = false;
      // Infection history to update
      new_infection_history = new_infection_history_mat(indiv,_);

      ///////////////////////////////////////////////////////
      // OPTION 1: Swap contents of a year for an individual
      ///////////////////////////////////////////////////////
      // If swap step      
      if(swap_step_option){
	//swap_proposals[indiv] += 1;
	loc1 = locs[j]; // Choose a location from age_mask to strain_mask
	loc2 = loc1 + floor(R::runif(-swap_distance,swap_distance)); // Perturb +/- swap_distance
      
	// If we have gone too far left or right, reflect at the boundaries
	while(loc2 < 0) loc2 += n_samp_length;
	if(loc2 >= n_samp_length) loc2 -= floor(loc2/n_samp_length)*n_samp_length;

	// Get onto right scale (starting at age mask)
	loc1 += age_mask[indiv] - 1;
	loc2 += age_mask[indiv] - 1;
	
	loc1_val_old = new_infection_history(loc1);
	loc2_val_old = new_infection_history(loc2);
	
	// Only proceed if we've actually made a change
	// If prior version 4, then prior doesn't change by swapping
	if(loc1_val_old != loc2_val_old){
	  lik_changed = true;
	  
	  if(!prior_on_total){
	    // Number of infections in that group in that time
	    m_1_old = n_infections(group_id,loc1);      
	    m_2_old = n_infections(group_id,loc2);
	  
	    // Swap contents
	    new_infection_history(loc1) = new_infection_history(loc2);
	    new_infection_history(loc2) = loc1_val_old;
	  
	    // Number alive is number alive overall in that time and group
	    n_1 = n_alive(group_id, loc1);
	    n_2 = n_alive(group_id, loc2);
	    
	    // Prior for new state
	    m_1_new = m_1_old - loc1_val_old + loc2_val_old;
	    m_2_new = m_2_old - loc2_val_old + loc1_val_old;

	    // Pre-compute these? 
	    prior_1_old = R::lbeta(m_1_old + alpha, n_1 - m_1_old + beta)-lbeta_const;
	    prior_2_old = R::lbeta(m_2_old + alpha, n_2 - m_2_old + beta)-lbeta_const;
	    prior_old = prior_1_old + prior_2_old;
	    
	    prior_1_new = R::lbeta(m_1_new + alpha, n_1 - m_1_new + beta)-lbeta_const;
	    prior_2_new = R::lbeta(m_2_new + alpha, n_2 - m_2_new + beta)-lbeta_const;
	    prior_new = prior_1_new + prior_2_new;
	  } else {
	    // Prior version 4
	    prior_old = prior_new = 0;
	  }
	}
	///////////////////////////////////////////////////////
	// OPTION 2: Add/remove infection
	///////////////////////////////////////////////////////
      } else {
	//add_proposals[indiv] += 1;
	year = locs[j] + age_mask[indiv] - 1;
	old_entry = new_infection_history(year);
	//Rcpp::Rcout << "Old entry: " << old_entry << std::endl;
	if(!prior_on_total){	
	  // Get number of individuals that were alive and/or infected in that year,
	  // less the current individual
	  // Number of infections in this year, less infection status of this individual in this year
	  m = n_infections(group_id, year) - old_entry;
	  n = n_alive(group_id, year) - 1;
	} else {
	  m = n_infected_group(group_id) - old_entry;
	  n = total_alive(group_id) - 1;
	}
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
	if(new_entry != old_entry){
	  lik_changed = true;
	}
      }
      ////////////////////////
      // If a change was made to the infection history,
      // calculate likelihood of new Z
      ////////////////////////
      if(solve_likelihood && lik_changed){
	// Calculate likelihood!
	indices = new_infection_history > 0;
	infection_times = circulation_times[indices];
	//if(infection_times.size() > 0){
	  infection_strain_indices_tmp = circulation_times_indices[indices];	  
	  // ====================================================== //
	  // =============== CHOOSE MODEL TO SOLVE =============== //
	  // ====================================================== //
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
	  } else if(alternative_wane_func){
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
	  } else {
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
	  //}
	if(use_titre_shifts){
	  add_measurement_shifts(predicted_titres, titre_shifts, 
				 start_index_in_data, end_index_in_data);
	}
	  
	// Now have all predicted titres for this individual calculated
	// Need to calculate likelihood of these titres... 
	new_prob = 0;
	
	// Go from first row in the data for this individual to up to the next one, accumlating
	// likelihood for this individual
	// For unique data

	proposal_likelihood_func(new_prob, predicted_titres, indiv, data, repeat_data, repeat_indices,
				 cum_nrows_per_individual_in_data, cum_nrows_per_individual_in_repeat_data,
				 log_const, den, max_titre);

      } else {
	old_prob = new_prob = old_probs[indiv];
      }
     
      //////////////////////////////
      // METROPOLIS-HASTINGS STEP
      //////////////////////////////
      if(swap_step_option){ 
	log_prob = std::min<double>(0.0, (new_prob+prior_new) - (old_prob+prior_old));
      } else {
	log_prob = std::min<double>(0.0, new_prob - old_prob);
      }
      rand1 = R::runif(0,1);
      if(lik_changed && log(rand1) < log_prob/temp){
	// Update the entry in the new matrix Z
	old_prob = new_prob;
	old_probs[indiv] = new_prob;

	// Carry out the swap
	if(swap_step_option){
	  tmp = new_infection_history_mat(indiv,loc1);
	  new_infection_history_mat(indiv,loc1) = new_infection_history_mat(indiv,loc2);
	  new_infection_history_mat(indiv,loc2) = tmp;
	  
	  // Update number of infections in the two swapped times
	  if(!prior_on_total){
	    n_infections(group_id, loc1) = m_1_new;
	    n_infections(group_id, loc2) = m_2_new;
	  }
	  // Don't need to update group infections if prior_on_total, as infections
	  // only move within an individual (so number in group stays same)
	} else {
	  new_infection_history_mat(indiv,year) = new_entry;	
	  // Update total number of infections in group/time
	  if(!prior_on_total){
	    n_infections(group_id, year) -= old_entry;
	    n_infections(group_id, year) += new_entry;
	  } else {
	    n_infected_group(group_id) = n_infected_group(group_id) - old_entry + new_entry;
	  }
	}
      }
    }
  }
  List ret;
  //Rcpp::Rcout << "Old probs at end: " << old_probs << std::endl;
  ret["old_probs"] = old_probs;
  ret["new_infection_history"] = new_infection_history_mat;
  return(ret);
}

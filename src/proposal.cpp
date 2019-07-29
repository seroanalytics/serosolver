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
				 const NumericVector& rand_ns,
				 const double& swap_propn) {
  // Copy input matrix
  arma::mat new_infection_history_mat = infection_history_mat;
  IntegerVector locs; // Locations to be updated
  arma::uvec locs1;
  arma::mat x;
  arma::mat y;
  IntegerVector samps;
  IntegerVector subset_samps;
  LogicalVector tmp_indices;
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
  int x_n_cols;
  
  double rand1;
  double ratio;
  
  // For each sampled individual
  for(int i = 0; i < sampled_indivs.size(); ++i){

    // Isolate that individual's infection histories
    indiv = sampled_indivs[i]-1;
    n_inf = n_infs[indiv];
    x = new_infection_history_mat.submat(indiv, age_mask[indiv]-1, indiv, strain_mask[indiv]-1);
    x_n_cols = x.n_cols;
    samps = seq(0,x_n_cols-1);
    // With some probability, add/remove infections or swap infections
    if(rand_ns[i] > swap_propn){
      n_samp_max = std::min(n_inf, x_n_cols); 
      // Sample N random locations
      locs = RcppArmadillo::sample(samps, n_samp_max, FALSE, NumericVector::create());
      locs1 = as<arma::uvec>(locs);
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
      tmp_indices = as<LogicalVector>(wrap(x));
      subset_samps = samps[tmp_indices];
      if(subset_samps.size() > 0){
	id1 = subset_samps(floor(R::runif(0,1)*subset_samps.size()));
	max_i_indiv = x.size();
	move_max = move_sizes[indiv];
	move = floor(R::runif(0,1)*2*move_max) - move_max;
	id2 = id1 + move;
	while(id2 < 0) id2 += max_i_indiv;
	while(id2 >= max_i_indiv) id2 -= max_i_indiv;
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
				   const NumericVector &antigenic_distances,
				   const NumericVector &data,
				   const NumericVector &repeat_data,
				   const IntegerVector &repeat_indices,
				   const NumericVector &titre_shifts,
				   IntegerVector proposal_iter,
				   IntegerVector accepted_iter,
				   IntegerVector proposal_swap,
				   IntegerVector accepted_swap,
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

  //Repeat data?
  bool repeat_data_exist = repeat_indices[0] >= 0;
  
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

 // Back boosting model
  double nu_long_recall;
  double nu_short_recall;
  double max_interference;
  double interference_gradient;
  double affinity_maturation;
  int back_boosting_type = theta["back_boosting"];
  bool back_boosting = back_boosting_type == 1;
  if(back_boosting){
    nu_long_recall = theta["nu_long"];
    nu_short_recall = theta["nu_short"];
    max_interference = theta["max_interference"];
    interference_gradient = theta["interference_gradient"];
    affinity_maturation = theta["affinity_maturation"];
  }
  
  // 3. If not using one of the specific mechanism functions, set the base_function flag to TRUE
  bool base_function = !(alternative_wane_func ||
			 titre_dependent_boosting ||
			 strain_dep_boost ||
			 back_boosting);
  
  // 4. Extra titre shifts
  bool use_titre_shifts = false;
  if(titre_shifts.size() == n_titres_total) use_titre_shifts = true;
  // ########################################################################
  ////Rcpp::Rcout << "Number strains: " << number_strains << std::endl;
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
    n_samp_length  = strain_mask[indiv] - age_mask[indiv] + 1; // How many times maximum can we sample from?
    ////Rcpp::Rcout << "N samp length: " << n_samp_length << std::endl;
    // If swap step, only doing one proposal for this individual
    if(swap_step_option){
      n_samp_max = 1;
      // Get this individual's infection history
      new_infection_history = new_infection_history_mat(indiv,_);

      // For the swap step, start by generating a vector from 0 to max     
      samps = seq(0, number_strains-1);    // Create vector across all potential infection times
      // Subset by entries that have an infection in them
      samps = samps[new_infection_history == 1];
      // Check if there were any infections. If not, we will skip
      if(samps.size() > 0){
	// If infections, start the samps vector such that 0 = age_mask
	samps = samps - age_mask[indiv] + 1;
	locs = RcppArmadillo::sample(samps, n_samp_max, FALSE, NumericVector::create());
      }
    } else {
      // Sample n_samp_length. Ths will be used to pull years from sample_years
      n_samp_max = std::min(n_years_samp, n_samp_length); // Use the smaller of these two numbers
      samps = seq(0, n_samp_length-1);    // Create vector from 0:length of alive years
      locs = RcppArmadillo::sample(samps, n_samp_max, FALSE, NumericVector::create());
    }
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
	prior_old = prior_new = 0;
	if(samps.size() > 0){	  
	  proposal_swap[indiv] += 1;
	  //swap_proposals[indiv] += 1;
	  loc1 = locs[j]; // Choose a location from age_mask to strain_mask
	  loc2 = loc1 + floor(R::runif(-swap_distance,swap_distance));

	  // If we have gone too far left or right, reflect at the boundaries
	  while(loc2 < 0){
	    // If gone negative, then reflect to the other side.
	    // ie. -1 becomes the last entry, -2 becomes the second last entry etc.
	    loc2 += n_samp_length;
	  }
	  while(loc2 >= n_samp_length){
	    loc2 -= n_samp_length;
	  }
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
	      ////Rcpp::Rcout << "Not prior on total" << std::endl;
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
	} 
	///////////////////////////////////////////////////////
	// OPTION 2: Add/remove infection
	///////////////////////////////////////////////////////
      } else {
	
	proposal_iter[indiv] += 1;
	//add_proposals[indiv] += 1;
	year = locs[j] + age_mask[indiv] - 1;
	old_entry = new_infection_history(year);
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
					  antigenic_map_long,
					  false);	  
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
					      false);	
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
						      false);
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
					   antigenic_map_long,
					   false);
	} else if(back_boosting) {
	  	titre_model_backboost_cpp(predicted_titres,
				  mu, mu_short,
				  wane, tau,
				  affinity_maturation,
				  nu_long_recall, nu_short_recall,
				  max_interference, interference_gradient,
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
				  antigenic_distances,
				  false);	  
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
					  antigenic_map_long,
					  false);
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
				 log_const, den, max_titre, repeat_data_exist);

      } else {
	old_prob = new_prob = old_probs[indiv];
      }
     
      //////////////////////////////
      // METROPOLIS-HASTINGS STEP
      //////////////////////////////
      if(swap_step_option){ 
	log_prob = std::min<double>(0.0, (new_prob+prior_new) - (old_prob+prior_old));
	////Rcpp::Rcout << "Log prob: " << log_prob << std::endl;
      } else {
	log_prob = std::min<double>(0.0, new_prob - old_prob);
      }
      rand1 = R::runif(0,1);
      if(lik_changed && log(rand1) < log_prob/temp){
	// Update the entry in the new matrix Z1
	old_prob = new_prob;
	old_probs[indiv] = new_prob;

	// Carry out the swap
	if(swap_step_option){
	  accepted_swap[indiv] += 1;
	  ////Rcpp::Rcout << "Carrying out swap" << std::endl;
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
	  accepted_iter[indiv] += 1;
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
  ////Rcpp::Rcout << "Old probs at end: " << old_probs << std::endl;
  ret["old_probs"] = old_probs;
  ret["new_infection_history"] = new_infection_history_mat;
  ret["proposal_iter"] = proposal_iter;
  ret["accepted_iter"] = accepted_iter;
  ret["proposal_swap"] = proposal_swap;
  ret["accepted_swap"] = accepted_swap;
  return(ret);
}



// [[Rcpp::export]]
double titre_data_fast_individual_base_indiv(const double &mu,
                                           const double &mu_short,
                                           const double &wane,
                                           const double &tau,
                                           const NumericVector &infection_times,
                                           const IntegerVector &infection_strain_indices_tmp,
                                           const int &measurement_strain_index,
                                           const double &sampling_time,
                                           const int &number_strains,
                                           const NumericVector &antigenic_map_short,
                                           const NumericVector &antigenic_map_long,
                                           bool boost_before_infection = false
){
  double predicted_titre=0;
  double time;
  double n_inf = 1.0;
  double wane_amount;
  double seniority;
  
  int n_titres;
  int max_infections = infection_times.size();
  int inf_map_index;
  int index;
  
  n_inf = 1.0;
  
  // Sum all infections that would contribute towards observed titres at this time
  for(int x = 0; x < max_infections; ++x){
    // Only go further if this sample happened after the infection
    if((boost_before_infection && sampling_time > infection_times[x]) ||
       (!boost_before_infection && sampling_time >= infection_times[x])){
      time = sampling_time - infection_times[x]; // Time between sample and infection
      wane_amount= MAX(0, 1.0 - (wane*time)); // Basic waning function
      seniority = MAX(0, 1.0 - tau*(n_inf - 1.0)); // Antigenic seniority
      inf_map_index = infection_strain_indices_tmp[x]; // Index of this infecting strain in antigenic map
      index = measurement_strain_index*number_strains + inf_map_index;
      // Find contribution to each measured titre from this infection
      predicted_titre += seniority * ((mu*antigenic_map_long[index]) +
        (mu_short*antigenic_map_short[index])*wane_amount);
      ++n_inf;
    }
  }
  return predicted_titre;
}

//[[Rcpp::export]]
double prob_x_given_z_titre_protection(const double &titre_p, const double &alpha1, const double &beta1,
			    const int &x, const int &z){
  if(x == 0 && z == 0){
    return 0;
  } else if(x == 0 && z == 1){
    return log(titre_p);
  } else if(x == 1 && z == 1){
    return log(1.0-titre_p);
  } else {
    return 0;
  }
  return 0;
}

//[[Rcpp::export]]
double prob_x_given_z_titre(const double &titre, const double &alpha1, const double &beta1,
			    const int &x, const int &z){
  double tmp_prob;
  if(x == 0 && z == 0){
    return 0;
  } else if(x == 0 && z == 1){
    tmp_prob = titre_protection_cpp(titre, alpha1, beta1);
    return log(tmp_prob);
  } else if(x == 1 && z == 1){
    tmp_prob = titre_protection_cpp(titre, alpha1, beta1);
    return log(1.0-tmp_prob);
  } else {
    return 0;
  }
  return 0;
}




//[[Rcpp::export]]
NumericMatrix calc_titre_inf_likelihoods(const NumericMatrix &titres,
				  const double &alpha1, const double &beta1,
				  const NumericMatrix &Xs, const NumericMatrix &Zs){
  int max_row = titres.nrow();
  int max_col = titres.ncol();
  NumericMatrix inf_liks(max_row, max_col);
  for(int i = 0; i < max_row; ++i){
    for(int j = 0; j < max_col; ++j){
      inf_liks(i,j) = prob_x_given_z_titre(titres(i,j), alpha1, beta1, Xs(i, j), Zs(i,j));
    }
  }
  return inf_liks;  
}
				  
//[[Rcpp::export]]
NumericMatrix calc_titre_inf_prob(const NumericMatrix &titres,
				  const double &alpha1, const double &beta1,
				  const NumericMatrix &Xs, const NumericMatrix &Zs){
  int max_row = titres.nrow();
  int max_col = titres.ncol();
  NumericMatrix inf_probs(max_row, max_col);
  for(int i = 0; i < max_row; ++i){
    for(int j = 0; j < max_col; ++j){
      inf_probs(i,j) = 1.0 - titre_protection_cpp(titres(i,j), alpha1, beta1);
    }
  }
  return inf_probs;  
}
	

// [[Rcpp::export]]
List inf_hist_prop_prior_immunity(const NumericVector &theta, // Model parameters
				  const IntegerMatrix &exposure_history_mat,  // Current exposure history
				  const IntegerMatrix &infection_history_mat, // Current infection history
				  const NumericMatrix &titre_prob_inf, // Probability of infection given exposure at each possible exposure time
				  const NumericVector &old_probs_1,
				  const IntegerVector &sampled_indivs,
				  const IntegerVector &n_years_samp_vec,
				  const IntegerVector &age_mask, // Age mask
				  const IntegerVector &strain_mask, // Age mask
				  const IntegerMatrix &n_alive, // No. of individuals alive each year/group
				  IntegerMatrix &n_exposures, // No. of infections in each year/group
				  const double &swap_propn,
				  const int &swap_distance,
				  const double &alpha_Z, // Alpha for prior
				  const double &beta_Z, // Beta for prior
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
				  IntegerVector proposal_iter,
				  IntegerVector accepted_iter,
				  IntegerVector proposal_swap,
				  IntegerVector accepted_swap,
				  const double temp=1,
				  bool solve_likelihood=true				   
				  ){
  // ########################################################################
  // Parameters to control indexing of data
  IntegerMatrix new_infection_history_mat(infection_history_mat);
  IntegerMatrix new_exposure_history_mat(exposure_history_mat);
  NumericMatrix new_titre_prob_inf = clone(titre_prob_inf);

  double min_titre = theta["min_titre"];
  int n_titres_total = data.size(); // How many titres are there in total?
  NumericVector predicted_titres(n_titres_total, min_titre); // Vector to store predicted titres
  double tmp_predicted_titre;
  NumericVector old_probs = clone(old_probs_1); // Create a copy of the current old probs
  
  // Variables related to solving likelihood and model as little as possible
  bool swap_step_option = true;
  bool lik_changed = true;
  //Repeat data?
  bool repeat_data_exist = repeat_indices[0] >= 0;
  
  // These quantities can be pre-computed
  int n_indivs = infection_history_mat.nrow();  // How many individuals are there in total?
  int number_strains = infection_history_mat.ncol(); // How many possible years are we interested in?
  int n_sampled = sampled_indivs.size(); // How many individuals are we actually investigating?
  int n_exposed = 0; // How many individuals were exposed?
    
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
  IntegerVector new_exposure_history(number_strains); // New proposed exposure history
  NumericVector new_titre_prob_inf_indiv(number_strains); // New proposed exposure history
  
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
  int next_year;

  
  int new_entry_z;
  int new_entry_x;
  int old_entry_z;
  int old_entry_x;
  
  int loc1, loc2, tmp; // Which indices are we looking at?
  int n_years_samp; // How many years to sample for this individual?

  int loc1_exposure_old, loc2_exposure_old;
  int loc1_infection_old, loc2_infection_old;
  double titre_prob_old, titre_prob_new;
  // ########################################################################

  // ########################################################################
  double m; // number of exposures in a given year
  double n; // number alive in a particular year

  double m_1_new, m_1_old,m_2_new,m_2_old;
  double n_1, n_2;
  double prior_1_old, prior_2_old, prior_1_new, prior_2_new, prior_new, prior_old;

  double prior_new_to_remove;
  double prior_old_to_remove;
  
  double rand1; // Store a random number
  double rand2;
  double ratio_exposure; // Store the gibbs ratio for 0 or 1 proposal
  double ratio_infection; // Store the gibbs ratio for 0 or 1 proposal

  double old_prob; // Likelihood of old number
  double new_prob; // Likelihood of new number
  double log_prob; // Likelihood ratio

  double lbeta_const = R::lbeta(alpha_Z, beta_Z);

  // For likelihood
  const double sd = theta["error"];
  const double den = sd*M_SQRT2;
  const double max_titre = theta["MAX_TITRE"];
  const double log_const = log(0.5);
  const double alpha_titre = theta["alpha_titre"];
  const double beta_titre = theta["beta_titre"];
  
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
  
  // Extra titre shifts
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
    //  Rcpp::Rcout << "Indiv: " << indiv << std::endl;
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

    // Get this individual's infection and exposure history and
    // titre probability terms
    new_exposure_history = new_exposure_history_mat(indiv,_);
    new_infection_history = new_infection_history_mat(indiv,_);
    new_titre_prob_inf_indiv = new_titre_prob_inf(indiv,_);

    // Assume that proposal hasn't changed likelihood until shown otherwise
    lik_changed = false;
    ///////////////////////////////////////////////////////
    // OPTION 1: Swap contents of a year for an individual
    ///////////////////////////////////////////////////////
    // If swap step      
    if(swap_step_option){
      samps = seq(0, number_strains-1);    // Create vector across all potential infection times
      indices = new_exposure_history > 0;
      samps = samps[indices];   // Subset by those indices that had infections

      // Old prior for infection state given exposure state
      prior_old = prior_new = 0;
      titre_prob_old = titre_prob_new = sum(new_titre_prob_inf_indiv);
      
      // Only proceed with swap step if there was a viable infection time to be chosen
      if(samps.size() > 0){
	proposal_swap[indiv] += 1;
	// If so, then we want to start "samps" from 
	samps = samps - age_mask[indiv] + 1;
	locs = RcppArmadillo::sample(samps, 1, FALSE, NumericVector::create());
	
	loc1 = locs[0]; // Choose a location from age_mask to strain_mask
	loc2 = loc1 + floor(R::runif(-swap_distance,swap_distance)); // Perturb +/- swap_distance

	// If we have gone too far left or right, reflect at the boundaries
	while(loc2 < 0){
	  // If gone negative, then reflect to the other side.
	  // ie. -1 becomes the last entry, -2 becomes the second last entry etc.
	  loc2 += n_samp_length;
	}
	while(loc2 >= n_samp_length){
	  loc2 -= n_samp_length;
	}
	
	
	// Get onto right scale (starting at age mask)
	loc1 += age_mask[indiv] - 1;
	loc2 += age_mask[indiv] - 1;

	// Get old exposure states at these two locations
	loc1_exposure_old = new_exposure_history(loc1);
	loc2_exposure_old = new_exposure_history(loc2);

	// Get old infection states at these locations
	loc1_infection_old = new_infection_history(loc1);
	loc2_infection_old = new_infection_history(loc2);

	// Only proceed if we've actually made a change
	//if((loc1_exposure_old != loc2_exposure_old) ||
	//   (loc1_infection_old != loc2_infection_old)) {
	  // Number of exposure in that group in that time
	  m_1_old = n_exposures(group_id,loc1);      
	  m_2_old = n_exposures(group_id,loc2);
	  
	  // Swap contents exposures
	  new_exposure_history(loc1) = new_exposure_history(loc2);
	  new_exposure_history(loc2) = loc1_exposure_old;

	  new_infection_history(loc1) = new_infection_history(loc2);
	  new_infection_history(loc2) = loc1_infection_old;

	  // Calculate titres at these infection times with new
	  // infection histories
	  // ***************************************************
	  indices = new_infection_history > 0;
	  infection_times = circulation_times[indices];
	  infection_strain_indices_tmp = circulation_times_indices[indices];
	  titre_prob_new = 0;
	  new_titre_prob_inf_indiv.fill(0);
	  for(int year_i = age_mask[indiv]-1; year_i < strain_mask[indiv]; ++year_i){
	    tmp_predicted_titre = titre_data_fast_individual_base_indiv(mu, mu_short,
									wane, tau,
									infection_times,
									infection_strain_indices_tmp,
									circulation_times_indices[year_i],
									circulation_times[year_i],
									number_strains,
									antigenic_map_short,
									antigenic_map_long,
									true);
	    // Calculate this given new_infection_history
	    old_entry_z = new_exposure_history[year_i];
	    old_entry_x = new_infection_history[year_i];
	    ratio_infection = titre_protection_cpp(tmp_predicted_titre, alpha_titre, beta_titre);
	    new_titre_prob_inf_indiv(year_i) = prob_x_given_z_titre_protection(ratio_infection,
									       alpha_titre, beta_titre,
									       old_entry_x, old_entry_z);
	  }
	  titre_prob_new = sum(new_titre_prob_inf_indiv);
	  // ***************************************************
	  // Number alive is number alive overall in that time and group
	  n_1 = n_alive(group_id, loc1);
	  n_2 = n_alive(group_id, loc2);
	    
	  // Prior for new state
	  m_1_new = m_1_old - loc1_exposure_old + loc2_exposure_old;
	  m_2_new = m_2_old - loc2_exposure_old + loc1_exposure_old;

	  // Old priors for exposure states
	  prior_1_old = R::lbeta(m_1_old + alpha_Z, n_1 - m_1_old + beta_Z)-lbeta_const;
	  prior_2_old = R::lbeta(m_2_old + alpha_Z, n_2 - m_2_old + beta_Z)-lbeta_const;
		
	  // Sum with titre probs
	  prior_old = prior_1_old + prior_2_old + titre_prob_old;

	  // New priors for exposure states
	  prior_1_new = R::lbeta(m_1_new + alpha_Z, n_1 - m_1_new + beta_Z)-lbeta_const;
	  prior_2_new = R::lbeta(m_2_new + alpha_Z, n_2 - m_2_new + beta_Z)-lbeta_const;

	  // Sum with titre probs
	  prior_new = prior_1_new + prior_2_new + titre_prob_new;
	  lik_changed = true;
	  //}
      }
      ///////////////////////////////////////////////////////
      // OPTION 2: Add/remove infection
      ///////////////////////////////////////////////////////
    } else {
      // For the swap step, start by generating a vector from 0 to max     
      n_samp_max = std::min(n_years_samp, n_samp_length); // Use the smaller of these two numbers
      samps = seq(0, n_samp_length);    // Create vector from 0:length of alive years
      locs = RcppArmadillo::sample(samps, n_samp_max, FALSE, NumericVector::create());
      std::sort(locs.begin(), locs.end());
      //    Rcpp::Rcout << "New_titre_prob_inf_indiv: " << new_titre_prob_inf_indiv << std::endl;
      //    Rcpp::Rcout << "Old exposure history: " << new_exposure_history << std::endl;
      //    Rcpp::Rcout << "Old infection history: " << new_infection_history << std::endl;
      //    Rcpp::Rcout << "Locs chosen: " << locs << std::endl;
      proposal_iter[indiv] += 1;
      
      // Old prior for infection state given exposure state
      prior_old = sum(new_titre_prob_inf_indiv);
      prior_new_to_remove = 0;
      prior_old_to_remove = 0;
      
      for(int j=0; j < n_samp_max; ++j){
	//	Rcpp::Rcout << "Loc selected: " << locs[j] << std::endl;
	indices = new_infection_history > 0;
	infection_times = circulation_times[indices];
	infection_strain_indices_tmp = circulation_times_indices[indices];     
	year = locs[j] + age_mask[indiv] - 1;

	//	Rcpp::Rcout << "Year selected: " << year << std::endl;;
	if(j < n_samp_max - 1){
	  next_year = locs[j+1]  + age_mask[indiv] - 1;
	} else {
	  next_year = strain_mask[indiv];
	}
	// Calculate infection prob likelihood for all times from the proposal
	// time onwards
	//for(int year_i = year; year_i < strain_mask[indiv]; ++year_i){
	// Location 1 titre and infection prob
	tmp_predicted_titre = titre_data_fast_individual_base_indiv(mu, mu_short,
								    wane, tau,
								    infection_times,
								    infection_strain_indices_tmp,
								    circulation_times_indices[year],
								    circulation_times[year],
								    number_strains,
								    antigenic_map_short,
								    antigenic_map_long,
								    true);
	// Calculate the probability of infection given exposure at this time
	ratio_infection = titre_protection_cpp(tmp_predicted_titre, alpha_titre, beta_titre);

	old_entry_z = new_exposure_history[year];
	old_entry_x = new_infection_history[year];

	//	Rcpp::Rcout << "Old entry z: " << old_entry_z << std::endl;
	//	Rcpp::Rcout << "Old entry x: " << old_entry_x << std::endl;
	
	/* ================= PROPOSE NEW ENTRY =================== */
	// Get number of individuals that were alive and/or infected in that year,
	// less the current individual
	// Number of infections in this year, less infection status of this individual in this year
	m = n_exposures(group_id, year) - old_entry_z;
	n = n_alive(group_id, year) - 1;
   
	// Work out proposal ratio - prior from alpha, beta and number of other infections
	ratio_exposure = (m + alpha_Z)/(n + alpha_Z + beta_Z);
	//	Rcpp::Rcout << "Ratio exposure: " << ratio_exposure << std::endl;
	//	Rcpp::Rcout << "Ratio infection: " << 1-ratio_infection << std::endl;
	
	// Propose 1 or 0 based on this ratio
	rand1 = R::runif(0,1);	
	if(rand1 < ratio_exposure){
	  new_entry_z = 1;
	  new_exposure_history(year) = 1;
	  rand2 = R::runif(0,1);
	  // Propose a new infection state given the exposure state.
	  // Just need prob infection at this time before the infection happens
	  if(rand2 < (1.0 - ratio_infection)){
	    new_entry_x = 1;
	    new_infection_history(year) = 1;
	  } else {
	    new_entry_x = 0;
	    new_infection_history(year) = 0;
	  }
	 
	} else {
	  new_entry_z = 0;
	  new_entry_x = 0;
	  new_exposure_history(year) = 0;
	  new_infection_history(year) = 0;
	}
	
	//Rcpp::Rcout << "New entry z: " << new_entry_z << std::endl;
	//Rcpp::Rcout << "New entry x: " << new_entry_x << std::endl;
	// Calculate probability of these new outcomes
	// Sampling from prior, so remove old prior prob for this state
	prior_old_to_remove += new_titre_prob_inf_indiv(year);
	new_titre_prob_inf_indiv(year) = prob_x_given_z_titre_protection(ratio_infection,
									   alpha_titre, beta_titre,
									   new_entry_x, new_entry_z);
	// Sampling from prior, so remove new prior prob for this state
	prior_new_to_remove += new_titre_prob_inf_indiv(year);
	
	year++;

	indices = new_infection_history > 0;
	infection_times = circulation_times[indices];
	infection_strain_indices_tmp = circulation_times_indices[indices];     
	// Update probabilities up to next location or end
	//for(int year_i = year; year_i < strain_mask[indiv]; ++year_i){
	while(year < next_year){

	  // Rcpp::Rcout << "Year: " << year << std::endl;
	  tmp_predicted_titre = titre_data_fast_individual_base_indiv(mu, mu_short,
								      wane, tau,
								      infection_times,
								      infection_strain_indices_tmp,
								      circulation_times_indices[year],
								      circulation_times[year],
								      number_strains,
								      antigenic_map_short,
								      antigenic_map_long,
								      true);
	  old_entry_z = new_exposure_history[year];
	  old_entry_x = new_infection_history[year];
	  // Calculate the probability of infection given exposure at this time
	  ratio_infection = titre_protection_cpp(tmp_predicted_titre, alpha_titre, beta_titre);
	  // Update prior prob for new state
	  
	  //Rcpp::Rcout << "Old titre prob inf: " << new_titre_prob_inf_indiv(year) << std::endl;
	  new_titre_prob_inf_indiv(year) = prob_x_given_z_titre_protection(ratio_infection,
									   alpha_titre, beta_titre,
									   old_entry_x, old_entry_z);
	  // Rcpp::Rcout << "Old entry z: " << old_entry_z << std::endl;
	  //  Rcpp::Rcout << "Old entry x: " << old_entry_x << std::endl;
	  //  Rcpp::Rcout << "New titre prob inf: " << new_titre_prob_inf_indiv(year) << std::endl;
	  year++;
	}
      }
      // Prior part of the acceptance ratio for the old state is the old prior, minus the reverse proposal
      // ratio
      prior_old = prior_old - prior_old_to_remove;
      prior_new = sum(new_titre_prob_inf_indiv) - prior_new_to_remove;
      //   Rcpp::Rcout << "Old prior: " << prior_old << std::endl;
      //   Rcpp::Rcout << "New prior: " << prior_new << std::endl;
      lik_changed = true;
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

      // Data likelihood
      // ====================================================== //
      // =============== CHOOSE MODEL TO SOLVE =============== //
      // ====================================================== //
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
				      antigenic_map_long,
				      false);
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
			       log_const, den, max_titre, repeat_data_exist);
      
      // Don't need to add titre mediated immunity probs, as sampling from these
    } else {
      old_prob = new_prob = old_probs[indiv];
    }
     
    //////////////////////////////
    // METROPOLIS-HASTINGS STEP
    //////////////////////////////
    if(swap_step_option){ 
      log_prob = std::min<double>(0.0, (new_prob+prior_new) - (old_prob+prior_old));
    } else {
      log_prob = std::min<double>(0.0, (new_prob+prior_new) - (old_prob+prior_old));
    }
    rand1 = R::runif(0,1);
    if(lik_changed && log(rand1) < log_prob/temp){
      // Update the entry in the new matrix Z
      old_prob = new_prob;
      old_probs[indiv] = new_prob;

      // Carry out the swap
      if(swap_step_option){
	accepted_swap[indiv] += 1;
	// Swap exposures
	// Old prior for infection state given exposure state
	new_exposure_history_mat(indiv,_) = new_exposure_history;
	new_infection_history_mat(indiv,_) = new_infection_history;
	new_titre_prob_inf(indiv,_) = new_titre_prob_inf_indiv;
	
	n_exposures(group_id, loc1) = m_1_new;
	n_exposures(group_id, loc2) = m_2_new;
	/*
	  tmp = new_exposure_history_mat(indiv,loc1);
	  new_exposure_history_mat(indiv,loc1) = new_exposure_history_mat(indiv,loc2);
	  new_exposure_history_mat(indiv,loc2) = tmp;

	  // Swap infections
	  tmp = new_infection_history_mat(indiv,loc1);
	  new_infection_history_mat(indiv,loc1) = new_infection_history_mat(indiv,loc2);
	  new_infection_history_mat(indiv,loc2) = tmp;
	*/		
      } else {
	accepted_iter[indiv] += 1;
	// Update total number of infections in group/time
	n_exposures(group_id, _) = n_exposures(group_id, _) - new_exposure_history_mat(indiv,_);	
	new_exposure_history_mat(indiv,_) = new_exposure_history;
	new_infection_history_mat(indiv,_) = new_infection_history;
	new_titre_prob_inf(indiv,_) = new_titre_prob_inf_indiv;	
	n_exposures(group_id, _) = n_exposures(group_id, _) + new_exposure_history_mat(indiv,_);
      }
    }
  }
  List ret;
  ret["titre_infection_probs"] = new_titre_prob_inf;
  ret["old_probs"] = old_probs;
  ret["new_infection_history"] = new_infection_history_mat;
  ret["new_exposure_history"] = new_exposure_history_mat;
  ret["proposal_iter"] = proposal_iter;
  ret["accepted_iter"] = accepted_iter;
  ret["proposal_swap"] = proposal_swap;
  ret["accepted_swap"] = accepted_swap;
  return(ret);
}

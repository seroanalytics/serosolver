#include <RcppArmadilloExtensions/sample.h>
#include "antibody_models_individual.h"
#include "likelihood_funcs.h"
#include "helpers.h"

//' Infection history proposal function
//' 
//' Proposes a new matrix of infection histories using a beta binomial proposal distribution. This particular implementation allows for n_infs epoch times to be changed with each function call. Furthermore, the size of the swap step is specified for each individual by proposal_inf_hist_distances.
//' @param infection_history_mat and RcppArmadillo matrix of infection histories, where rows represent individuals and columns represent potential infection times. The contents should be a set of 1s (presence of infection) and 0s (absence of infection)
//' @param sampled_indivs IntegerVector, indices of which individuals to resample. Note that this is indexed from 1 (ie. as if passing straight from R)
//' @param age_mask IntegerVector, for each individual gives the first column in the infection history matrix that an individual could have been exposed to indexed from 1. ie. if alive for the whole period, entry would be 1. If alive for the 11th epoch, entry would be 11.
//' @param strain_mask IntegerVector, for each individual gives the last column in the infection history matrix that an individual could have been exposed to indexed from 1. ie. if their last serum sample was in the 40th epoch, entry would be 40
//' @param proposal_inf_hist_distances IntegerVector, how far can a swap step sample from specified for each individual
//' @param n_infs IntegerVector, how many infections to add/remove/swap with each proposal step for each individual
//' @param shape1 double, shape1 (alpha) parameter of the beta binomial
//' @param shape2 double, shape2 (beta) parameter of the beta binomial
//' @param rand_ns NumericVector, a vector of random numbers for each sampled individual. The idea is to pre-specify whether an individual experiences an add/remove step or a swap step to avoid random number sampling in C++
//' @return a matrix of 1s and 0s corresponding to the infection histories for all individuals
//' @export
//' @family infection_history_proposal
// [[Rcpp::export]]
arma::mat inf_hist_prop_prior_v3(arma::mat infection_history_mat, 
				 const IntegerVector& sampled_indivs, 
				 const IntegerVector& age_mask,
				 const IntegerVector& sample_mask,
				 const IntegerVector& proposal_inf_hist_distances, 
				 const IntegerVector& n_infs,
				 double shape1, 
				 double shape2, 
				 const NumericVector& rand_ns,
				 const double& proposal_inf_hist_indiv_swap_ratio) {
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
    x = new_infection_history_mat.submat(indiv, age_mask[indiv]-1, indiv, sample_mask[indiv]-1);
    x_n_cols = x.n_cols;
    samps = seq(0,x_n_cols-1);
    // With some probability, add/remove infections or swap infections
    if(rand_ns[i] > proposal_inf_hist_indiv_swap_ratio){
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
        ratio = (shape1 + k)/(shape1 + shape2 + n);
        rand1 = R::runif(0,1);
	// With probability 'ratio', add a 1. ie. if many 1s already, less likely
	// to add more 1s depending on shape1 and shape2
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
	move_max = proposal_inf_hist_distances[indiv];
	move = floor(R::runif(0,1)*2*move_max) - move_max;
	id2 = id1 + move;
	while(id2 < 0) id2 += max_i_indiv;
	while(id2 >= max_i_indiv) id2 -= max_i_indiv;
	tmp = x[id1];
	x[id1] = x[id2];
	x[id2] = tmp;
      }
    }
    new_infection_history_mat.submat(indiv, age_mask[indiv]-1, indiv,  sample_mask[indiv]-1) = x;
  }
  return(new_infection_history_mat);
}




//' Infection history gibbs proposal
//'
//' Generates a new infection history matrix and corresponding individual likelihoods, using a gibbs sampler from the infection history prior. See \code{\link{inf_hist_prop_prior_v3}}, as inputs are very similar.
//' @param theta NumericVector, the named model parameters used to solve the model
//' @param infection_history_mat IntegerMatrix the matrix of 1s and 0s corresponding to individual infection histories
//' @param likelihoods_pre_proposal NumericVector, the current likelihoods for each individual before proposing new infection histories
//' @param sampled_indivs IntegerVector, indices of sampled individuals
//' @param n_times_samp_vec int, for each individual, how many time periods to resample infections for?
//' @param age_mask IntegerVector, length of the number of individuals, with indices specifying first time period that an individual can be infected (indexed from 1, such that a value of 1 allows an individual to be infected in any time period)
//' @param sample_mask IntegerVector, length of the number of individuals, with indices specifying last time period that an individual can be infected (ie. last time a sample was taken)
//' @param n_alive IntegerMatrix, number of columns is the number of time periods that an individual could be infected, giving the number of individual alive in each time period. Number of rows is the number of distinct groups.
//' @param n_infections IntegerMatrix, the number of infections in each year (columns) for each group (rows)
//' @param n_infected_group IntegerVector, the total number of infections across all times in each group
//' @param prior_lookup arma::cube, the pre-computed lookup table for the beta prior on infection histories, dimensions are number of infections, time, and group
//' @param proposal_inf_hist_indiv_swap_ratio double, gives the proportion of proposals that will be swap steps (ie. swap contents of two cells in infection_history rather than adding/removing infections)
//' @param swap_distance int, in a swap step, how many time steps either side of the chosen time period to swap with
//' @param shape1 double, shape1 (alpha) parameter for beta prior on infection probability
//' @param shape2 double, shape2 (beta) parameter for beta prior on infection probability
//' @param possible_exposure_times NumericVector, the times that individuals could be infected
//' @param possible_exposure_times_indices IntegerVector, indexing vector from 0:(number of exposure times-1)
//' @param sample_times NumericVector, the vector of real times that samples were taken
//' @param rows_per_indiv_in_samples IntegerVector, How many rows in antibody data correspond to each individual, sample and repeat?
//' @param cum_nrows_per_individual_in_data IntegerVector, How many rows in the antibody data correspond to each individual?
//' @param cum_nrows_per_individual_in_repeat_data IntegerVector, For the repeat data (ie. already calculated these antibody levels), how many rows in the antibody data correspond to each individual?
//' @param nrows_per_blood_sample IntegerVector, Split the sample times and runs for each individual
//' @param group_id_vec IntegerVector, vector with 1 entry per individual, giving the group ID of that individual
//' @param biomarker_id_indices IntegerVector, For each antibody measurement, corresponding entry in antigenic map
//' @param antigenic_map_long arma::mat, the collapsed cross reactivity map for long term boosting, after multiplying by sigma1, see \code{\link{create_cross_reactivity_vector}}
//' @param antigenic_map_short arma::mat, the collapsed cross reactivity map for short term boosting, after multiplying by sigma2, see \code{\link{create_cross_reactivity_vector}}
//' @param antigenic_distances NumericVector matching the dimensions of antigenic_map_long and antigenic_map_short, but with the raw antigenic distances between strains
//' @param antibody_data NumericVector, data for all individuals for the first instance of each calculated antibody level
//' @param antibody_data_repeats NumericVector, the repeat antibody data for all individuals (ie. do not solve the same antibody level twice)
//' @param repeat_indices IntegerVector, which index in the main data vector does each entry in repeat_data correspond to ie. which calculated antibody level in predicted_antibody_levels should be used for each observation?
//' @param measurement_shifts NumericVector, if length matches the length of \code{data}, adds these as measurement shifts to the antibody levels. If lengths do not match, is not used.
//' @param proposal_iter IntegerVector, vector with entry for each individual, storing the number of infection history add/remove proposals for each individual.
//' @param accepted_iter IntegerVector, vector with entry for each individual, storing the number of accepted infection history add/remove proposals for each individual.
//' @param proposal_swap IntegerVector, vector with entry for each individual, storing the number of proposed infection history swaps
//' @param accepted_swap IntegerVector, vector with entry for each individual, storing the number of accepted infection history swaps
//' @param total_alive IntegerVector, giving the total number of potential infection events for each group. This only applies to prior version 4. If set to a vector of values -1, then this is ignored.
//' @param temp double, temperature for parallel tempering MCMC
//' @param solve_likelihood bool, if FALSE does not solve likelihood when calculating acceptance probability
//' @param data_type int, defaults to 1 for discretized, bounded data. Set to 2 for continuous, bounded data
//' @return an R list with 6 entries: 1) the vector replacing likelihoods_pre_proposal, corresponding to the new likelihoods per individual; 2) the matrix of 1s and 0s corresponding to the new infection histories for all individuals; 3-6) the updated entries for proposal_iter, accepted_iter, proposal_swap and accepted_swap.
//' @export
//' @family infection_history_proposal
// [[Rcpp::export]]
List inf_hist_prop_prior_v2_and_v4(
        
        
        const NumericVector &theta, //All model parameters
        const IntegerVector &unique_theta_indices, //Indices for each model parameter type
        const IntegerVector &unique_biomarker_groups, // Vector of unique observation types
        
        
				   const IntegerMatrix &infection_history_mat,  // Current infection history
				   const IntegerVector &infection_history_mat_indices,
				   const NumericVector &likelihoods_pre_proposal,
				   const IntegerVector &sampled_indivs,
				   const IntegerVector &n_times_samp_vec,
				   const IntegerVector &age_mask, // Age mask
				   const IntegerVector &sample_mask, // Age mask
				   const IntegerMatrix &n_alive, // No. of individuals alive each year/group
				   IntegerMatrix &n_infections, // No. of infections in each year/group
				   IntegerVector &n_infected_group,
				   const arma::cube &prior_lookup,
				   const double &proposal_inf_hist_indiv_swap_ratio,
				   const int &swap_distance,
				   const bool &propose_from_prior,
				   const double &shape1, // shape1 for prior
				   const double &shape2, // shape2 for prior
				   const NumericVector &possible_exposure_times,
				   const IntegerVector &possible_exposure_times_indices,
				   const NumericVector &sample_times,
				   
				   const IntegerVector &type_data_start, // For each individual, which entry in the unique(antibody _data[,c("individual","biomarker_group)]) data frame is their first?
				   const IntegerVector &biomarker_groups, // For each unique sample/individual index, what observation type is it?
				   
				   const IntegerVector &sample_data_start, 
				   const IntegerVector &antibody_data_start, 
				   const IntegerVector &nrows_per_sample, // Split the sample times for each individual
				   
				   const IntegerVector &cum_nrows_per_individual_in_data, // How many rows in the antibody data correspond to each individual? By biomarker_group
				   const IntegerVector &cum_nrows_per_individual_in_repeat_data, // How many rows in the repeat antibody data correspond to each individual? By biomarker_group
				   
				   const IntegerVector &group_id_vec, // Which group does each individual belong to?
				   const IntegerVector &biomarker_id_indices, // For each measurement, corresponding entry in antigenic map
				   const IntegerVector &start_level_indices, // For each individual/biomarker group/biomarker id combo, we need to know the starting antibody level
				   const NumericVector &starting_antibody_levels,
				   const NumericVector &births,
				   
				   const arma::mat &antigenic_map_long, // Now a matrix of antigenic maps
				   const arma::mat &antigenic_map_short,
				   
				   const NumericVector &antigenic_distances,
				   
				   const NumericVector &antibody_data,
				   const NumericVector &antibody_data_repeats,
				   
				   const int &n_measurements_total,
				   
				   const IntegerVector &repeat_indices, // Which indices in the main antibody data do we use for repeat measurements?
				   const bool &repeat_data_exist,
				   
				   const NumericVector &measurement_shifts,
				   IntegerVector proposal_iter,
				   IntegerVector accepted_iter,
				   IntegerVector proposal_swap,
				   IntegerVector accepted_swap,
				   IntegerMatrix overall_swap_proposals,
				   IntegerMatrix overall_add_proposals,
				   
				   const NumericVector time_sample_probs,
				   const IntegerVector &total_alive,
				   
				   const IntegerVector &data_types,
				   const NumericVector &obs_weights,
				   const double temp=1,
				   bool solve_likelihood=true){
  // ########################################################################
  // Parameters to control indexing of data
  IntegerMatrix new_infection_history_mat(clone(infection_history_mat)); // Can this be avoided? Create a copy of the inf hist matrix
  
  List ret; // To return

  NumericVector predicted_antibody_levels(n_measurements_total); // Vector to store antibody levels
  NumericVector likelihoods_pre_proposal_tmp = clone(likelihoods_pre_proposal); // Create a copy of the current old probs
  
  // Variables related to solving likelihood and model as little as possible
  bool swap_step_option = true;
  bool lik_changed = false;
  
  // These quantities can be pre-computed
  int number_possible_exposures = possible_exposure_times.size(); // infection_history_mat.ncol(); // How many possible years are we interested in?
  int n_sampled = sampled_indivs.size(); // How many individuals are we actually investigating?
  
  // Using prior version 2 or 4?
  bool prior_on_total = total_alive(0) > 0;

  // To track how far through the larger vectors we move for each individual
  int biomarker_group=0;
  double obs_weight=1.0;
  int data_type = 1;
  int type_start;
  int type_end;
  int start_index_in_samples;
  int end_index_in_samples;
  int start_index_in_data;
  int end_index_in_data;

  int group_id; // Vector of group IDs for each individual
 
  IntegerVector new_infection_history(number_possible_exposures); // New proposed infection history
  IntegerVector infection_history(number_possible_exposures); // Old infection history
  LogicalVector indices;
  NumericVector infection_times; // Tmp store infection times for this infection history, combined with indices
  IntegerVector infection_times_indices_tmp; // Tmp store which index in antigenic map these infection times relate to
  IntegerVector use_indices;
  
   
  // ########################################################################
  // Parameters related to infection history sampling
  int indiv; // Index of the individual under consideration
  IntegerVector samps; // Variable vector to sample from
  IntegerVector samps_shifted;
  IntegerVector locs; // Vector of locations that were sampled
  // As each individual has a different number of years to sample from,
  // need to extract relative proportions and re-weight
  NumericVector tmp_loc_sample_probs;
  int year; // Index of year being updated
  int n_samp_max; // Maximum number of years to sample for individual
  int n_samp_length; // Number of years that COULD be sampled for this individual

  int new_entry = 0;
  int old_entry = 0;
  int loc1 = 0;
  int loc2 = 0;
  int tmp = 0; // Which indices are we looking at?
  int n_times_samp; // How many years to sample for this individual?
  int loc1_val_old, loc2_val_old;
  // ########################################################################

  // ########################################################################
  double m; // number of infections in a given year
  double n; // number alive in a particular year
  //double n_1; // Number alive in time period 1 
  //double n_2; // Number alive in time period 2
  
  double m_1_new, m_1_old,m_2_new,m_2_old;
  double prior_1_old, prior_2_old, prior_1_new,prior_2_new,prior_new,prior_old;

  double rand1; // Store a random number
  double ratio; // Store the gibbs ratio for 0 or 1 proposal

  double old_prob; // Likelihood of old number
  double new_prob; // Likelihood of new number
  double log_prob; // Likelihood ratio

  //double lbeta_const = R::lbeta(shape1, shape2);

  // ====================================================== //
  // =============== SETUP MODEL PARAMETERS =============== //
  // ====================================================== //
  // 1. Extract general parameters that apply to all models
  // Pull out model parameters so only need to allocate once
  int n_types = unique_biomarker_groups.size();
  int n_theta = unique_theta_indices.size();
  
  // For likelihood functions
  NumericVector sds(n_types);
  NumericVector dens(n_types);
  const double log_const = log(0.5);
  NumericVector den2s(n_types);
  
  // Base model parameters
  NumericVector boost_long_parameters(n_types);
  NumericVector boost_short_parameters(n_types);
  NumericVector boost_delay_parameters(n_types);
  NumericVector wane_short_parameters(n_types);
  NumericVector wane_long_parameters(n_types);
  NumericVector antigenic_seniority_parameters(n_types);
  
  NumericVector max_measurements(n_types);
  NumericVector min_measurements(n_types);

  int boost_long_index = unique_theta_indices("boost_long");
  int boost_short_index = unique_theta_indices("boost_short");
  int boost_delay_index = unique_theta_indices("boost_delay");
  int wane_short_index = unique_theta_indices("wane_short");
  int wane_long_index = unique_theta_indices("wane_long");
  int antigenic_seniority_index = unique_theta_indices("antigenic_seniority");
  int error_index = unique_theta_indices("obs_sd");
  
  int min_index = unique_theta_indices("min_measurement");
  int max_index = unique_theta_indices("max_measurement");
  
  // Antibody-dependent boosting function
  IntegerVector antibody_dependent_boosting(n_types);
  NumericVector gradients(n_types); 
  NumericVector boost_limits(n_types);   
  
  int antibody_dependent_boosting_index = unique_theta_indices("antibody_dependent_boosting");
  int gradient_index = -1;
  int boost_limit_index = -1;
  
  // Create vectors of model parameters for each of the observation types
  for(int x = 0; x < n_types; ++x){
    boost_long_parameters(x) = theta(boost_long_index + x*n_theta);
    boost_short_parameters(x) = theta(boost_short_index + x*n_theta);
    boost_delay_parameters(x) = theta(boost_delay_index + x*n_theta);
    wane_short_parameters(x) = theta(wane_short_index + x*n_theta);
    wane_long_parameters(x) = theta(wane_long_index + x*n_theta);
    antigenic_seniority_parameters(x) = theta(antigenic_seniority_index + x*n_theta);
      min_measurements(x) = theta(min_index + x*n_theta);
      max_measurements(x) = theta(max_index + x*n_theta);
      
      // For likelihood functions
      sds(x) = theta(error_index + x*n_theta);
      dens(x) = sds(x)*M_SQRT2; // Constant for the discretized normal distribution
      den2s(x) = log(sds(x)*2.50662827463); // Constant for the normal distribution
      
      // Titre dependent boosting
      antibody_dependent_boosting(x) = theta(antibody_dependent_boosting_index+ x*n_theta);
      //Rcpp::Rcout << "Titre dependent: " << antibody_dependent_boosting(x) << std::endl;
      if(antibody_dependent_boosting(x) == 1) {
          gradient_index = unique_theta_indices("gradient");
          boost_limit_index = unique_theta_indices("boost_limit");  
          gradients(x) = theta(gradient_index + x*n_theta);
          boost_limits(x) = theta(boost_limit_index + x*n_theta);
      }
  }

  // 4. Extra titre shifts
  bool use_measurement_shifts = false;
  if(measurement_shifts.size() == n_measurements_total) use_measurement_shifts = true;
  
  // ########################################################################
  // For each individual
  for(int i = 0; i < n_sampled; ++i){
    // Which proposal step to take and do we need to calculate the likelihood    
    swap_step_option = R::runif(0,1) < proposal_inf_hist_indiv_swap_ratio;
    
    // Get index, group and current likelihood of individual under consideration
    indiv = sampled_indivs(i)-1;
    //Rcpp::Rcout << "indiv: " << indiv+1 << std::endl << std::endl;
    
    group_id = group_id_vec(indiv);
    old_prob = likelihoods_pre_proposal(indiv);
    //Rcpp::Rcout << "Old prob: " << old_prob << std::endl << std::endl;
    
    // Time sampling control
    n_times_samp = n_times_samp_vec(indiv); // How many times are we intending to resample for this individual?
    n_samp_length  = sample_mask(indiv) - age_mask(indiv) + 1; // How many times maximum can we sample from?
    
    // If swap step, only doing one proposal for this individual
    if(swap_step_option){
      n_samp_max = 1;
      // Get this individual's infection history
      new_infection_history = new_infection_history_mat(indiv,_);
    } else {
      // Sample n_samp_length. Ths will be used to pull years from sample_years
      n_samp_max = std::min(n_times_samp, n_samp_length); // Use the smaller of these two numbers
    }
    samps = seq(0, n_samp_length-1);    // Create vector from 0:length of alive years

    // Extract time sampling probabilities and re-normalise
    samps_shifted = samps + age_mask(indiv) - 1;
   
    tmp_loc_sample_probs = time_sample_probs[samps_shifted];
    // Re-normalise
    tmp_loc_sample_probs = tmp_loc_sample_probs/sum(tmp_loc_sample_probs);
    locs = RcppArmadillo::sample(samps, n_samp_max, false, tmp_loc_sample_probs);

    
    // For each selected infection history entry
    for(int j = 0; j < n_samp_max; ++j){
        //Rcpp::Rcout << "Updating infection history entry: " << j << std::endl;
      // Assume that proposal hasn't changed likelihood until shown otherwise
      lik_changed = false;
      // Infection history to update
      new_infection_history = new_infection_history_mat(indiv,_);
      //Rcpp::Rcout << "Infection history before: " << new_infection_history<< std::endl;
      
      ///////////////////////////////////////////////////////
      // OPTION 1: Swap contents of a year for an individual
      ///////////////////////////////////////////////////////
      // If swap step
      prior_old = prior_new = 0;
      if(swap_step_option){
    	loc1 = locs(j); // Choose a location from age_mask to sample_mask
    	loc2 = loc1 + floor(R::runif(-swap_distance,swap_distance+1));
    
    	// If we have gone too far left or right, reflect at the boundaries
    	
    	while(loc2 < 0){
    	  // If gone negative, then reflect to the other side.
    	  // ie. -1 becomes the last entry, -2 becomes the second last entry etc.
    	  loc2 += n_samp_length;
    	}
    	while(loc2 >= n_samp_length){
    	  loc2 -= n_samp_length;
    	}
	
    	// Try bounce rather than reflect to other side
    	//if(loc2 < 0) loc2 = -loc2;
    	//if(loc2 >= n_samp_length) loc2 = n_samp_length - loc2 + n_samp_length - 2;
    	
    	// Get onto right scale (starting at age mask)
    	loc1 += age_mask(indiv) - 1;
    	loc2 += age_mask(indiv) - 1;
    	  
    	loc1_val_old = new_infection_history(loc1);
    	loc2_val_old = new_infection_history(loc2);
    
    	overall_swap_proposals(indiv,loc1)++;
    	overall_swap_proposals(indiv,loc2)++;
	
    	// Only proceed if we've actually made a change
    	// If prior version 4, then prior doesn't change by swapping
    	if(loc1_val_old != loc2_val_old){
    	  lik_changed = true;
    	  proposal_swap(indiv) += 1;
    	  if(!prior_on_total){
    	    // Number of infections in that group in that time
    	      m_1_old = n_infections(group_id,loc1);      
    	      m_2_old = n_infections(group_id,loc2);
    	  
    	      // Swap contents
    	      new_infection_history(loc1) = new_infection_history(loc2);
    	      new_infection_history(loc2) = loc1_val_old;
    	  
    	      // Number alive is number alive overall in that time and group
    	      //n_1 = n_alive(group_id, loc1);
    	      //n_2 = n_alive(group_id, loc2);
    	    
    	      // Prior for new state
    	      m_1_new = m_1_old - loc1_val_old + loc2_val_old;
    	      m_2_new = m_2_old - loc2_val_old + loc1_val_old;
    
    	      
    	      prior_1_old = prior_lookup(m_1_old, loc1, group_id);
    	      prior_2_old = prior_lookup(m_2_old, loc2, group_id);
    	      prior_old = prior_1_old + prior_2_old;
    	      
    	      prior_1_new = prior_lookup(m_1_new, loc1, group_id);
    	      prior_2_new = prior_lookup(m_2_new, loc2, group_id);
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
    	year = locs(j) + age_mask(indiv) - 1;
    	old_entry = new_infection_history(year);
    	overall_add_proposals(indiv,year)++;
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

    	if(propose_from_prior){
            // Is there room to add another infection?
    	  // Work out proposal ratio - prior from shape1, shape2 and number of other infections
    	  ratio = (m + shape1)/(n + shape1 + shape2);
    
    	  // Propose 1 or 0 based on this ratio
    	  rand1 = R::runif(0,1);	

    	  if(rand1 < ratio){
    	    new_entry = 1;
    	    new_infection_history(year) = 1;
    	  } else {
    	    new_entry = 0;
    	    new_infection_history(year) = 0;
    	  }
    	} else {
    	  if(old_entry == 0) {
    	    new_entry = 1;
    	    new_infection_history(year) = 1;
    	  } else {
    	    new_entry = 0;
    	    new_infection_history(year) = 0;
    	  }
    	  m_1_old = m + old_entry;
    	  m_1_new = m + new_entry;
    	  
    	  prior_old = prior_lookup(m_1_old, year, group_id);
    	  prior_new = prior_lookup(m_1_new, year, group_id);
    	}
    	
    	if(new_entry != old_entry){
    	  lik_changed = true;
    	  proposal_iter(indiv) += 1;		
    	}
      }
      ////////////////////////
      // If a change was made to the infection history,
      // calculate likelihood of new Z
      ////////////////////////
      
      if(solve_likelihood && lik_changed){
        //if(TRUE){
         //Rcpp::Rcout << "Infection history after change: " << new_infection_history<< std::endl;
          
    	// Calculate likelihood!
    	indices = new_infection_history > 0;
      use_indices =infection_history_mat_indices[indices];
    	infection_times = possible_exposure_times[use_indices];
    	infection_times_indices_tmp = possible_exposure_times_indices[use_indices];	  

    	//Rcpp::Rcout << infection_times_indices_tmp << std::endl;
    	
    	// Start end end location of the type_data matrix
    	type_start = type_data_start(indiv);
    	type_end = type_data_start(indiv+1)-1;
    	
    	//Rcpp::Rcout << "Type start: " << type_start << std::endl;
    	//Rcpp::Rcout << "Type end: " << type_end << std::endl << std::endl;
    	
    	// ====================================================== //
    	// =============== CHOOSE MODEL TO SOLVE =============== //
    	// ====================================================== //
    	// For each observation type solved for this individual
	    new_prob = 0;

    	for(int index = type_start; index <= type_end; ++index){
    	    //Rcpp::Rcout << "index: " << index << std::endl;
    	    biomarker_group = biomarker_groups(index)-1;
    	    data_type = data_types(biomarker_group);
    	    obs_weight = obs_weights(biomarker_group);
    	 
    	    // Now have all antibody levels for this individual calculated
    	    // Need to calculate likelihood of these measurements... 
    	    start_index_in_samples = sample_data_start(index);
    	    end_index_in_samples = sample_data_start(index+1)-1;
    	    start_index_in_data = antibody_data_start(start_index_in_samples);
    	    end_index_in_data = antibody_data_start(end_index_in_samples+1)-1;
    	    
    	    
        	if (antibody_dependent_boosting(biomarker_group)) {
        	  antibody_dependent_boosting_model_individual(predicted_antibody_levels, 
                                               boost_long_parameters(biomarker_group), 
                                               boost_short_parameters(biomarker_group),
                                               wane_short_parameters(biomarker_group), 
                                               antigenic_seniority_parameters(biomarker_group),
                                               gradients(biomarker_group), 
                                               boost_limits(biomarker_group),
        					      infection_times,
        					      infection_times_indices_tmp,
        					      biomarker_id_indices,
        					      sample_times,
        					      start_index_in_samples,
        					      end_index_in_samples,
        					      start_index_in_data,
        					      nrows_per_sample,
        					      number_possible_exposures,
        					      antigenic_map_short.colptr(biomarker_group),
        					      antigenic_map_long.colptr(biomarker_group),
        					      false);	
        	} else {
        	  //Rcpp::Rcout << "Min measurement: " << min_measurements(biomarker_group) << std::endl;
        	  antibody_data_model_individual_new(
        	    predicted_antibody_levels, 
        	    starting_antibody_levels,
        	    births,
        	    boost_long_parameters(biomarker_group), 
        	    boost_short_parameters(biomarker_group),
        	    boost_delay_parameters(biomarker_group),
        	    wane_short_parameters(biomarker_group), 
        	    wane_long_parameters(biomarker_group), 
        	    antigenic_seniority_parameters(biomarker_group),
        	    infection_times,
        	    infection_times_indices_tmp,
        	    biomarker_id_indices,
        	    start_level_indices,
        	    sample_times,
        	    start_index_in_samples,
        	    end_index_in_samples,
        	    start_index_in_data,
        	    nrows_per_sample,
        	    number_possible_exposures,
        	    antigenic_map_short.colptr(biomarker_group),
        	    antigenic_map_long.colptr(biomarker_group),
        	    false,
        	    min_measurements(biomarker_group));
        	  /*
        	  antibody_data_model_individual(predicted_antibody_levels, 
                                          boost_long_parameters(biomarker_group), 
                                          boost_short_parameters(biomarker_group),
                                          boost_delay_parameters(biomarker_group),
                                          wane_short_parameters(biomarker_group), 
                                          wane_long_parameters(biomarker_group), 
                                           antigenic_seniority_parameters(biomarker_group),
        					  infection_times,
        					  infection_times_indices_tmp,
        					  biomarker_id_indices,
        					  sample_times,
        					  start_index_in_samples,
        					  end_index_in_samples,
        					  start_index_in_data,
        					  nrows_per_sample,
        					  number_possible_exposures,
        					  antigenic_map_short.colptr(biomarker_group),
        					  antigenic_map_long.colptr(biomarker_group),
        					  false);
        	   */
        	}
        	if(use_measurement_shifts){
        	  add_measurement_shifts(predicted_antibody_levels, measurement_shifts, 
        				 start_index_in_data, end_index_in_data);
        	}
        	// Go from first row in the data for this individual to up to the next one, accumulating
        	// likelihood for this individual
        	// For unique data
            // Data_type 2 is continuous, bounded data
              if(data_type==2){
                proposal_likelihood_func_continuous(new_prob, predicted_antibody_levels, 
                                                    index, 
                                                    antibody_data, 
                                                    antibody_data_repeats, 
                                                    repeat_indices,
                                                    cum_nrows_per_individual_in_data, 
                                                    cum_nrows_per_individual_in_repeat_data,
                                                    log_const, 
                                                    sds(biomarker_group), dens(biomarker_group), den2s(biomarker_group), 
                                                    max_measurements(biomarker_group), 
                                                    min_measurements(biomarker_group), 
                                                    repeat_data_exist,
                                                    obs_weight);
            	
              } else {
                // Data_type 1 is discretized, bounded data
                proposal_likelihood_func(new_prob, predicted_antibody_levels, 
                                         index, 
                                         antibody_data, 
                                         antibody_data_repeats, 
                                         repeat_indices,
                                         cum_nrows_per_individual_in_data, 
                                         cum_nrows_per_individual_in_repeat_data,
                                         log_const, 
                                         dens(biomarker_group), 
                                         max_measurements(biomarker_group), 
                                         min_measurements(biomarker_group), 
                                         repeat_data_exist,
                                         obs_weight);
              }
    	}
      } else {
        old_prob = new_prob = likelihoods_pre_proposal_tmp(indiv);
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
    	// Update the entry in the new matrix Z1
    	old_prob = new_prob;
    	likelihoods_pre_proposal_tmp(indiv) = new_prob;
    
    	// Carry out the swap
    	if(swap_step_option){
    	  accepted_swap(indiv) += 1;
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
    	    
    	    accepted_iter(indiv) += 1;
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
  ret["likelihoods_pre_proposal_tmp"] = likelihoods_pre_proposal_tmp;
  ret["new_infection_history"] = new_infection_history_mat;
  ret["proposal_iter"] = proposal_iter;
  ret["accepted_iter"] = accepted_iter;
  ret["proposal_swap"] = proposal_swap;
  ret["accepted_swap"] = accepted_swap;
  ret["overall_swap_proposals"] = overall_swap_proposals;
  ret["overall_add_proposals"] = overall_add_proposals;
  return(ret);
}

#include <RcppArmadilloExtensions/sample.h>
#include "infection_model.h"
#include "time.h"
#include "iostream"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

#ifndef MAX
#define MAX(a,b) ((a) < (b) ? (b) : (a)) // define MAX function for use later
#endif



// [[Rcpp::export]]
NumericVector subset_nullable_vector(const Nullable<NumericVector> &x, int index1, int index2) {
  if(x.isNotNull()){
    NumericVector y = as<NumericVector>(x)[Range(index1, index2)];
    return y;
  } else {
    NumericVector y(1);
    return y;
  }
}


//' Gibbs sampling of infection histories
//'
//' Proposes a new matrix of infection histories by sampling from the prior on an individual's infection presence/absence in a particular time period, conditional on all other individuals in that time period. This allows us to integrate out the infection probability term for each time period. Should look at \code{\link{create_post_func}} for more details about the input parameters.
//' @param pars NumericVector, the model parameters used to solve the model, extract likelihood function parameters and alpha/beta
//' @param infHist IntegerMatrix the matrix of 1s and 0s corresponding to individual infection histories
//' @param indivSampPropn double, what proportion of individuals to resample in this proposal step
//' @param n_years_samp int, for each individual, how many time periods to resample infections for?
//' @param ageMask IntegerVector, length of the number of individuals, with indices specifying first time period that an individual can be infected (indexed from 1, such that a value of 1 allows an individual to be infected in any time period)
//' @param n_alive IntegerVector, length of the number of time periods that an individual could be infected, giving the number of individual alive in each time period
//' @param swapPropn double, what proportion of proposals should be swap steps (ie. swap contents of two cells in infHist rather than adding/removing infections)
//' @param swapDistance int, in a swap step, how many time steps either side of the chosen time period to swap with
//' @param alpha double, alpha parameter for beta prior on infection probability
//' @param beta double, beta parameter for beta prior on infection probability
//' @param circulationTimes NumericVector, the times that each strain circulated
//' @param circulationMapIndices IntegerVector, indexing vector from 1:number of strains
//' @param samplingTimes NumericVector, the vector of real times that samples were taken
//' @param indicesTitreDataSamples IntegerVector, How many rows in titre data correspond to each individual, sample and repeat?
//' @param indicesTitreDataOverall IntegerVector, How many rows in the titre data correspond to each individual?
//' @param indicesSamples IntegerVector, Split the sample times and runs for each individual
//' @param measuredMapIndices IntegerVector, For each titre measurement, corresponding entry in antigenic map
//' @param antigenicMapLong NumericVector, the collapsed cross reactivity map for long term boosting, after multiplying by sigma1
//' @param antigenicMapShort NumericVector, the collapsed cross reactivity map for short term boosting, after multiplying by sigma2
//' @param data NumericVector, the titre data for all individuals
//' @param to_add Nullable<NumericVector>, optional vector of measurement shifts to apply to all titres
//' @param additional_arguments Nullable<List>, optional list to pass down to titre solving function
//' @param DOBs NumericVector, vector of ages for each individual
//' @param solve_likelihood bool, if FALSE does not solve likelihood when calculating acceptance probability
//' @param total_alive int, if a positive number, uses this rather than the vector of alive individuals
//' @param temp double, temperature for parallel tempering MCMC
//' @return a matrix of 1s and 0s corresponding to the infection histories for all individuals
// [[Rcpp::export]]
IntegerMatrix infection_history_proposal_gibbs(const NumericVector& pars, // Model parameters
					       const IntegerMatrix& infHist,  // Current infection history
					       double indivSampPropn, // Proportion of individuals to resample
					       const IntegerVector& n_years_samp_vec, // Number of years to resample for each year
 					       const IntegerVector& ageMask, // Age mask
 					       const IntegerVector& strainMask, // Age mask
					       const IntegerVector& n_alive, // Number of individuals alive in each year
					       double swapPropn,
					       int swapDistance,
					       double alpha, // Alpha for prior
					       double beta, // Beta for prior
					       const NumericVector &circulationTimes,
					       const IntegerVector &circulationMapIndices,
					       const NumericVector &samplingTimes,
					       const IntegerVector &indicesTitreDataSample, // How many rows in titre data correspond to each individual, sample and repeat?
					       const IntegerVector &indicesTitreDataOverall, // How many rows in the titre data correspond to each individual?
					       const IntegerVector &indicesSamples, // Split the sample times and runs for each individual
					       const IntegerVector &measuredMapIndices, // For each titre measurement, corresponding entry in antigenic map
					       const NumericVector &antigenicMapLong, 
					       const NumericVector &antigenicMapShort,
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
  IntegerMatrix newInfHist(infHist);
  int n_indivs = newInfHist.nrow();
  int n_strains = newInfHist.ncol();
  // These indices allow us to step through the titre data vector
  // as if it were a matrix ie. number of rows for each individual
  // at a time
  int startIndexSamples;
  int endIndexSamples;
  int startIndexData;
  int endIndexData;
  // ########################################################################

  // ########################################################################
  IntegerVector infs_in_year;
  IntegerVector samps;
  IntegerVector sample_years;
  IntegerVector locs;
  IntegerVector proposedIndivHist;
  IntegerVector indivHist;
  IntegerVector swapSize;
  // ########################################################################

  // ########################################################################
  int indiv; // Index of the individual under consideration
  int year; // Index of year being updated
  double DOB; // The DOB of the individual for solving the model, if needed
  int maxYears = newInfHist.ncol(); // Number of epochs of circulation
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
    n_infected = sum(infHist);
  }
  // ########################################################################
  // For each individual
  for(int i = 1; i <= n_indivs; ++i){
    // Indexing for data upkeep
    startIndexSamples = indicesSamples[i-1];
    endIndexSamples = indicesSamples[i] - 1;

    startIndexData = indicesTitreDataOverall[i-1];
    endIndexData = indicesTitreDataOverall[i] - 1;
    if (to_add.isNotNull()){
      to_add_tmp = subset_nullable_vector(to_add, startIndexData, endIndexData);
    }
    
    // Choose whether to sample this individual or not
    if(R::runif(0,1) < indivSampPropn){
      // Index of this individual
      indiv = i-1;
      DOB = DOBs[indiv];
      n_years_samp = n_years_samp_vec[indiv];
      // Make vector of year indices to sample from
      // These are the indices in the matrix Z
      sample_years = seq(ageMask[indiv]-1,strainMask[indiv]-1);
      // Sample the minimum of either the number of years alive, or 
      // the number of years that are asked to be changed
      n_samp_max = sample_years.size();
      samps = seq(0, n_samp_max-1);    // Create vector from 0:length of alive years
      n_samp_max = std::min(n_years_samp, n_samp_max);
      // Swap contents of a year for an individual
      if(R::runif(0,1) < swapPropn){
	proposedIndivHist = newInfHist(indiv,_);
	indivHist = newInfHist(indiv,_);
	  
	loc1 = samps(floor(R::runif(0,samps.size())));
	swapSize = seq(-swapDistance,swapDistance);
	loc2 = loc1 + swapSize(floor(R::runif(0,swapSize.size())));

	while(loc2 < 0) loc2 = loc2 + samps.size();
	if(loc2 >= samps.size()) loc2 = loc2 - floor(loc2/samps.size())*samps.size();

	loc1 = sample_years(loc1);
	loc2 = sample_years(loc2);

	if(newInfHist(indiv,loc1) != newInfHist(indiv, loc2)){
	  tmp = proposedIndivHist(loc1);
	  proposedIndivHist(loc1) = proposedIndivHist(loc2);
	  proposedIndivHist(loc2) = tmp;

	  if(!prior_on_total){
	    // Prior for previous state
	    m_1_old = sum(newInfHist(_,loc1));      
	    m_2_old = sum(newInfHist(_,loc2));
	    
	    // Number alive is number alive overall
	    n_1_old = n_alive(loc1);
	    n_2_old = n_alive(loc2);

	    prior_1_old = R::lbeta(m_1_old + alpha, n_1_old - m_1_old + beta)-R::lbeta(alpha,beta);
	    prior_2_old = R::lbeta(m_2_old + alpha, n_2_old - m_2_old + beta)-R::lbeta(alpha,beta);
	    prior_old = prior_1_old + prior_2_old;

	    // Prior for new state
	    m_1_new = sum(newInfHist(_,loc1)) - newInfHist(indiv,loc1) + proposedIndivHist(loc1);      
	    m_2_new = sum(newInfHist(_,loc2)) - newInfHist(indiv,loc2) + proposedIndivHist(loc2);

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
						  indivHist, 
						  circulationTimes, 
						  circulationMapIndices,
						  samplingTimes[Range(startIndexSamples, endIndexSamples)], 
						  indicesTitreDataSample[Range(startIndexSamples,endIndexSamples)], 
						  measuredMapIndices[Range(startIndexData,endIndexData)],
						  antigenicMapLong, 
						  antigenicMapShort,  
						  n_strains,
						  data[Range(startIndexData,endIndexData)],
						  to_add_tmp,
						  DOB,
						  additional_arguments);

	    new_prob = likelihood_data_individual(pars, proposedIndivHist,
						  circulationTimes, circulationMapIndices,
						  samplingTimes[Range(startIndexSamples, endIndexSamples)], 
						  indicesTitreDataSample[Range(startIndexSamples,endIndexSamples)], 
						  measuredMapIndices[Range(startIndexData,endIndexData)],
						  antigenicMapLong, 
						  antigenicMapShort,  
						  n_strains,
						  data[Range(startIndexData,endIndexData)],
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
	    tmp = newInfHist(indiv,loc1);
	    newInfHist(indiv,loc1) = newInfHist(indiv,loc2);
	    newInfHist(indiv,loc2) = tmp;
	  }
	}
      } else {
	// Sample n_samp_real years from 0:length. Ths will be used to pull years from
	// sample_years
	// Note sampling indices in the individual's infection history, not the matrix Z
	locs = RcppArmadillo::sample(samps, n_samp_max, FALSE, NumericVector::create());
	// For each sampled year
	/*int k = sum(newInfHist);
	  int nn = total_alive;
	  for(int ii = 0; ii < locs.size(); ++ii){
	  k -= newInfHist(indiv, sample_years(locs[ii]));
	  }
	  nn -= locs.size();
	  proposedIndivHist = newInfHist(indiv,_);
	  indivHist = newInfHist(indiv,_);
	*/
	for(int j = 0; j < n_samp_max; ++j){
	  proposedIndivHist = newInfHist(indiv,_);
	  indivHist = newInfHist(indiv,_);
	
	  //  ratio = (k + alpha)/(nn + alpha + beta);
	  // Which year under consideration
	  // Get index for the matrix Z (not index for individual's inf hist)
	  // Note that locs[j] is the index in the individual's inf hist
	  year = sample_years(locs[j]);
	  	  // alpha = alphas[year];
	  // beta = betas[year];

	  if(prior_on_total){
	    old_entry = newInfHist(indiv, year);
	    // NUmber infected overall
	    m = n_infected - old_entry;
	    // Number alive is number alive overall
	    n = total_alive - 1;
	  } else {
	    // Get number of individuals that were alive and/or infected in that year,
	    // less the current individual
	    // Number of infections in this year, less infection status of this individual in this year
	    m = sum(newInfHist(_,year)) - newInfHist(indiv,year);     
	    n = n_alive(year) - 1;
	  }
	  // Work out proposal ratio - prior from alpha, beta and number of other infections
	  ratio = (m + alpha)/(n + alpha + beta);
	  // Propose 1 or 0 based on this ratio
	  rand1 = R::runif(0,1);
	  //if(indivRatioRands[i] < ratio){
	  if(rand1 < ratio){
	    new_entry = 1;
	    //k++;
	    proposedIndivHist(year) = 1;
	  } else {
	    new_entry = 0;
	    proposedIndivHist(year) = 0;
	  }
	  //nn++;
	
	  // If proposing a change, need to check likelihood ratio
	  if(new_entry != newInfHist(indiv,year)){
	    if(solve_likelihood){
	      old_prob = likelihood_data_individual(pars, indivHist, circulationTimes, circulationMapIndices,
						    samplingTimes[Range(startIndexSamples, endIndexSamples)], 
						    indicesTitreDataSample[Range(startIndexSamples,endIndexSamples)], 
						    measuredMapIndices[Range(startIndexData,endIndexData)],
						    antigenicMapLong, 
						    antigenicMapShort,  
						    n_strains,
						    data[Range(startIndexData,endIndexData)],
						    to_add_tmp,
						    DOB,
						    additional_arguments);
	      new_prob = likelihood_data_individual(pars, proposedIndivHist,
						    circulationTimes, circulationMapIndices,
						    samplingTimes[Range(startIndexSamples, endIndexSamples)], 
						    indicesTitreDataSample[Range(startIndexSamples,endIndexSamples)], 
						    measuredMapIndices[Range(startIndexData,endIndexData)],
						    antigenicMapLong, 
						    antigenicMapShort,  
						    n_strains,
						    data[Range(startIndexData,endIndexData)],
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
	      newInfHist(indiv, year) = new_entry;
	      //newInfHist(indiv, _) = proposedIndivHist;
	    }
	  }
	}
      }
    }
  }
  return(newInfHist);
}

//' Fast infection history proposal function
//' 
//' Proposes a new matrix of infection histories using a beta binomial proposal distribution. This particular implementation allows for nInfs epoch times to be changed with each function call. Furthermore, the size of the swap step is specified for each individual by moveSizes.
//' @param infHist and RcppArmadillo matrix of infection histories, where rows represent individuals and columns represent potential infection times. The contents should be a set of 1s (presence of infection) and 0s (absence of infection)
//' @param sampledIndivs IntegerVector, indices of which individuals to resample. Note that this is indexed from 1 (ie. as if passing straight from R)
//' @param ageMask IntegerVector, for each individual gives the first column in the infection history matrix that an individual could have been exposed to indexed from 1. ie. if alive for the whole period, entry would be 1. If alive for the 11th epoch, entry would be 11.
//' @param moveSizes IntegerVector, how far can a swap step sample from specified for each individual
//' @param nInfs IntegetVector, how many infections to add/remove/swap with each proposal step for each individual
//' @param alpha double, alpha parameter of the beta binomial
//' @param beta double, beta parameter of the beta binomial
//' @param randNs NumericVector, a vector of random numbers for each sampled individual. The idea is to pre-specify whether an individual experiences an add/remove step or a swap step to avoid random number sampling in C++
//' @return a matrix of 1s and 0s corresponding to the infection histories for all individuals
// [[Rcpp::export]]
arma::mat inf_hist_prop_cpp(arma::mat infHist, 
			    const IntegerVector& sampledIndivs, 
			    const IntegerVector& ageMask,
			    const IntegerVector& strainMask,
			    const IntegerVector& moveSizes, 
			    const IntegerVector& nInfs,
			    double alpha, 
			    double beta, 
			    const NumericVector& randNs) {
  // Copy input matrix
  arma::mat newInfHist = infHist;
  IntegerVector locs; // Locations to be updated
  arma::uvec locs1;
  arma::mat x;
  arma::mat y;
  IntegerVector samps;
  int maxI_indiv;
  int indiv;
  int k;
  int nInf;
  int n;
  int moveMax;
  int move;
  int id1;
  int id2;
  int tmp;
  int n_samp_max;

  double rand1;
  double ratio;
  
  // For each sampled individual
  for(int i = 0; i < sampledIndivs.size(); ++i){

    // Isolate that individual's infection histories
    indiv = sampledIndivs[i]-1;
    nInf = nInfs[indiv];
    x = newInfHist.submat(indiv, ageMask[indiv]-1, indiv, strainMask[indiv]-1);
    samps = seq_len(x.n_cols);
    n_samp_max = samps.size();
    n_samp_max = std::min(nInf, n_samp_max);
    
    // With 50% probability, add/remove infections or swap infections
    if(randNs[i] < 1.0/2.0){
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
      maxI_indiv = x.size();
      IntegerVector moves;
      for(int j = 0; j < n_samp_max; j++){
        id1 = floor(R::runif(0,1)*x.size());
        moveMax = moveSizes[indiv];
        move = floor(R::runif(0,1)*2*moveMax) - moveMax;
        id2 = id1 + move;
	while(id2 < 0) id2 += maxI_indiv;
	if(id2 >= maxI_indiv) id2 = (id2 % maxI_indiv);
	tmp = x[id1];
        x[id1] = x[id2];
        x[id2] = tmp;
      }
    }
    newInfHist.submat(indiv, ageMask[indiv]-1, indiv,  strainMask[indiv]-1) = x;
  }
  return(newInfHist);
}



// @export
// [[Rcpp::export]]
List infection_history_proposal_gibbs_fast(const NumericVector &theta, // Model parameters
					   const IntegerMatrix &infection_history_mat,  // Current infection history
					   const NumericVector &old_probsA,
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
					   const double temp=1,
					   bool solve_likelihood=true
					   ){
  // ########################################################################
  // Parameters to control indexing of data
  IntegerMatrix new_infection_history_mat(infection_history_mat); // Can this be avoided? Create a copy of the inf hist matrix
  int n_titres_total = data.size(); // How many titres are there in total?
  NumericVector predicted_titres(n_titres_total); // Vector to store predicted titres
  NumericVector old_probs = clone(old_probsA); // Create a copy of the current old probs

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
  int max_infections; // Tmp store of number of infections for this individual
  int n_titres; // Tmp store of number of titres to predict for this sample
 
  double sampling_time; // Tmp store time that blood sample taken
  double time; // Tmp store time between sample and exposure

  IntegerVector new_infection_history(number_strains); // New proposed infection history
  IntegerVector infection_history(number_strains); // Old infection history
  LogicalVector indices;

  NumericVector infection_times; // Tmp store infection times for this infection history, combined with indices
  IntegerVector infection_strain_indices_tmp; // Tmp store which index in antigenic map these infection times relate to

  double mu = theta["mu"];
  double mu_short = theta["mu_short"];
  double wane = theta["wane"];
  double tau = theta["tau"];
  double wane_amount;
  double seniority;
  double n_inf;
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
  // ########################################################################
  
  // ########################################################################
  // For each individual
  for(int i = 0; i < n_sampled; ++i){
    // Get index of individual under consideration and their current likelihood
    indiv = sampled_indivs[i]-1;
    old_prob = old_probsA[indiv];  

    // Indexing for data upkeep
    index_in_samples = rows_per_indiv_in_samples[indiv];
    end_index_in_samples = rows_per_indiv_in_samples[indiv+1] - 1;
    number_samples = end_index_in_samples - index_in_samples;      

    n_years_samp = n_years_samp_vec[indiv]; // How many years are we intending to resample from for this individual?
    n_samp_length  = strain_mask[indiv] - age_mask[indiv]; // How many years maximum can we sample from?
    n_samp_max = std::min(n_years_samp, n_samp_length); // Use the smaller of these two numbers

    // Swap contents of a year for an individual
    if(R::runif(0,1) < swap_propn){
      // Indexing for data upkeep
      start_index_in_data = cum_nrows_per_individual_in_data[indiv];
      start_index_in_repeat_data = cum_nrows_per_individual_in_repeat_data[indiv];

      new_infection_history = new_infection_history_mat(indiv,_);

      loc1 = floor(R::runif(0,n_samp_length)); // Choose a location from age_mask to strain_mask
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
	  max_infections = infection_times.size();
	  
	  // Don't need to check that no. infections >0, as would only swap if there was a 0 and a 1
	  for(int j = index_in_samples; j <= end_index_in_samples; ++j){
	    sampling_time = sample_times[j];
	    n_inf = 1.0;	
	    n_titres = nrows_per_blood_sample[j];
	
	    end_index_in_data = start_index_in_data + n_titres;
	    tmp_titre_index = start_index_in_data;

	    // ==================================
	    // THE MODEL
	    // ==================================
	    for(int x = 0; x < max_infections; ++x){
	      if(sampling_time >= infection_times[x]){
		time = sampling_time - infection_times[x];
		wane_amount= MAX(0, 1.0 - (wane*time));
		seniority = MAX(0, 1.0 - tau*(n_inf - 1.0));
		inf_map_index = infection_strain_indices_tmp[x];

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
	  // Now have all predicted titres for this individual calculated
	  // Need to calculate likelihood of these titres... 
	  new_prob = 0;

	  // NEED SOMETHING TO DEAL WITH REPEATS
	  // Go from first row in the data for this individual to up to the next one, accumlating
	  // likelihood for this individual
	  // For unique data
	  //Rcpp::Rcout << "Starting index in titre data: " << cum_nrows_per_individual_in_data[indiv] << std::endl;
	  for(int x = cum_nrows_per_individual_in_data[indiv]; x < cum_nrows_per_individual_in_data[indiv+1]; ++x){
	    //Rcpp::Rcout << "Predicted titre " << x << ": " << predicted_titres[x] << std::endl;
	    if(data[x] <= max_titre && data[x] > 0.0){
	      new_prob += log_const + log((erf((data[x] + 1.0 - predicted_titres[x]) / den) -
					   erf((data[x]     - predicted_titres[x]) / den)));    
	    } else if(data[x] > max_titre) {
	      new_prob += log_const + log(1.0 + erf((max_titre - predicted_titres[x])/den));
	    } else {
	      new_prob += log_const + log(1.0 + erf((1.0 - predicted_titres[x])/den));
	    }
	  }

	  // =====================
	  // Do something for repeat data here
	  for(int x = cum_nrows_per_individual_in_repeat_data[indiv]; x < cum_nrows_per_individual_in_repeat_data[indiv+1]; ++x){
	    if(repeat_data[x] <= max_titre && repeat_data[x] > 0.0){
	      new_prob += log_const + log((erf((repeat_data[x] + 1.0 - predicted_titres[repeat_indices[x]]) / den) -
					   erf((repeat_data[x]     - predicted_titres[repeat_indices[x]]) / den)));    
	    } else if(repeat_data[x] > max_titre) {
	      new_prob += log_const + log(1.0 + erf((max_titre - predicted_titres[repeat_indices[x]])/den));
	    } else {
	      new_prob += log_const + log(1.0 + erf((1.0 - predicted_titres[repeat_indices[x]])/den));
	    }
	  }
	  // Need to erase the predicted titre data...
	  for(int x = cum_nrows_per_individual_in_data[indiv]; x < cum_nrows_per_individual_in_data[indiv+1]; ++x){
	    predicted_titres[x] = 0;
	  }
	  

	  // =====================
	} else {
	  old_prob = new_prob = 0;
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
	start_index_in_data = cum_nrows_per_individual_in_data[indiv];
	start_index_in_repeat_data = cum_nrows_per_individual_in_repeat_data[indiv];

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
	      max_infections = infection_times.size();	  
	      infection_strain_indices_tmp = circulation_times_indices[indices];

	      for(int q = index_in_samples; q <= end_index_in_samples; ++q){
		sampling_time = sample_times[q];
		n_inf = 1.0;	
		n_titres = nrows_per_blood_sample[q];
	
		end_index_in_data = start_index_in_data + n_titres;
		tmp_titre_index = start_index_in_data;

		// ==================================
		// THE MODEL
		// ==================================
		for(int x = 0; x < max_infections; ++x){
		  if(sampling_time >= infection_times[x]){
		    time = sampling_time - infection_times[x];
		    wane_amount= MAX(0, 1.0 - (wane*time));
		    seniority = MAX(0, 1.0 - tau*(n_inf - 1.0));
		    inf_map_index = infection_strain_indices_tmp[x];

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
	    // Now have all predicted titres for this individual calculated
	    // Need to calculate likelihood of these titres... 
	    new_prob = 0;

	    // NEED SOMETHING TO DEAL WITH REPEATS
	    // Go from first row in the data for this individual to up to the next one, accumlating
	    // likelihood for this individual
	    // For unique data
	    for(int x = cum_nrows_per_individual_in_data[indiv]; x < cum_nrows_per_individual_in_data[indiv+1]; ++x){
	      if(data[x] <= max_titre && data[x] > 0.0){
		new_prob += log_const + log((erf((data[x] + 1.0 - predicted_titres[x]) / den) -
					     erf((data[x]     - predicted_titres[x]) / den)));    
	      } else if(data[x] > max_titre) {
		new_prob += log_const + log(1.0 + erf((max_titre - predicted_titres[x])/den));
	      } else {
		new_prob += log_const + log(1.0 + erf((1.0 - predicted_titres[x])/den));
	      }
	    }
	    // =====================
	    // Do the same thing but for repeat data
	    for(int x = cum_nrows_per_individual_in_repeat_data[indiv]; x < cum_nrows_per_individual_in_repeat_data[indiv+1]; ++x){
	      if(repeat_data[x] <= max_titre && repeat_data[x] > 0.0){
		new_prob += log_const + log((erf((repeat_data[x] + 1.0 - predicted_titres[repeat_indices[x]]) / den) -
					     erf((repeat_data[x]     - predicted_titres[repeat_indices[x]]) / den)));    
	      } else if(repeat_data[x] > max_titre) {
		new_prob += log_const + log(1.0 + erf((max_titre - predicted_titres[repeat_indices[x]])/den));
	      } else {
		new_prob += log_const + log(1.0 + erf((1.0 - predicted_titres[repeat_indices[x]])/den));
	      }
	    }
	    // Need to erase the predicted titre data...
	    for(int x = cum_nrows_per_individual_in_data[indiv]; x < cum_nrows_per_individual_in_data[indiv+1]; ++x){
	      predicted_titres[x] = 0;
	    }	  
	    // =====================
	  } else {
	    old_prob = new_prob = 0;
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

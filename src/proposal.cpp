#include <RcppArmadilloExtensions/sample.h>
#include "infection_model.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

#define MAX(a,b) ((a) < (b) ? (b) : (a)) // define MAX function for use later

//' Fast infection history proposal function
//' @return a matrix of 1s and 0s corresponding to the infection histories for all individuals
// [[Rcpp::export]]
IntegerMatrix infection_history_proposal_gibbs(NumericVector pars, // Model parameters
					       IntegerMatrix infHist,  // Current infection history
					       double indivSampPropn, // Proportion of individuals to resample
					       int n_years_samp, // Number of years to resample for each year
					       IntegerVector ageMask, // Age mask
					       IntegerVector n_alive, // Number of individuals alive in each year
					       double swapPropn,
					       int swapDistance,
					       double alpha, // Alpha for prior
					       double beta, // Beta for prior
					       NumericVector circulationTimes,
					       IntegerVector circulationMapIndices,
					       NumericVector samplingTimes,
					       IntegerVector indicesTitreDataSample, // How many rows in titre data correspond to each individual, sample and repeat?
					       IntegerVector indicesTitreDataOverall, // How many rows in the titre data correspond to each individual?
					       IntegerVector indicesSamples, // Split the sample times and runs for each individual
					       IntegerVector measuredMapIndices, // For each titre measurement, corresponding entry in antigenic map
					       NumericVector antigenicMapLong, 
					       NumericVector antigenicMapShort,
					       NumericVector data
					       ){
  // ########################################################################
  // Parameters to control indexing of data
  IntegerMatrix newInfHist = infHist;
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
  // 1. Copy infection history matrix
 
  //Rcpp::Rcout << sum(newInfHist) << std::endl;
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
  int maxYears = newInfHist.ncol(); // Number of epochs of circulation
  int n_samp_max; // Maximum number of years to sample
  int new_entry;
  int loc1, loc2, tmp;
  // ########################################################################

  // ########################################################################
  double m; // number of infections in a given year
  double n; // number alive in a particular year

  double rand1; // Store a random number
  double ratio; // Store the gibbs ratio for 0 or 1 proposal

  double old_prob; // Likelihood of old number
  double new_prob; // Likelihood of new number
  double log_prob; // Likelihood ratio
  // ########################################################################

  // For each individual
  for(int i = 1; i <= n_indivs; ++i){
    // Indexing for data upkeep
    startIndexSamples = indicesSamples[i-1];
    endIndexSamples = indicesSamples[i] - 1;

    startIndexData = indicesTitreDataOverall[i-1];
    endIndexData = indicesTitreDataOverall[i] - 1;
  
    // Choose whether to sample this individual or not
    if(R::runif(0,1) < indivSampPropn){
      // Index of this individual
      indiv = i-1;   
      // Make vector of year indices to sample from
      // These are the indices in the matrix Z
      sample_years = seq(ageMask[indiv]-1,maxYears-1);
      // Sample the minimum of either the number of years alive, or 
      // the number of years that are asked to be changed
      n_samp_max = sample_years.size();
      samps = seq(0, n_samp_max-1);    // Create vector from 0:length of alive years
      n_samp_max = std::min(n_years_samp, n_samp_max);
   
      // Swap contents of a year for an individual
      if(R::runif(0,1) > swapPropn){
	proposedIndivHist = newInfHist(indiv,_);
	indivHist = newInfHist(indiv,_);
	  
	loc1 = samps(floor(R::runif(0,samps.size())));
	swapSize = seq(-swapDistance,swapDistance);
	loc2 = loc1 + swapSize(floor(R::runif(0,swapSize.size())));

	if(loc2 < 0) loc2 = loc2 + samps.size();
	if(loc2 >= samps.size()) loc2 = loc2 - floor(loc2/samps.size())*samps.size();

	if(newInfHist(indiv,loc1) != newInfHist(indiv, loc2)){
	  tmp = proposedIndivHist(loc1);
	  proposedIndivHist(loc1) = proposedIndivHist(loc2);
	  proposedIndivHist(loc2) = tmp;
	  old_prob = likelihood_data_individual(pars, indivHist, circulationTimes, circulationMapIndices,
						samplingTimes[Range(startIndexSamples, endIndexSamples)], 
						indicesTitreDataSample[Range(startIndexSamples,endIndexSamples)], 
						measuredMapIndices[Range(startIndexData,endIndexData)],
						antigenicMapLong, 
						antigenicMapShort,  
						n_strains,
						data[Range(startIndexData,endIndexData)]);
	  new_prob = likelihood_data_individual(pars, proposedIndivHist,
						circulationTimes, circulationMapIndices,
						samplingTimes[Range(startIndexSamples, endIndexSamples)], 
						indicesTitreDataSample[Range(startIndexSamples,endIndexSamples)], 
						measuredMapIndices[Range(startIndexData,endIndexData)],
						antigenicMapLong, 
						antigenicMapShort,  
						n_strains,
						data[Range(startIndexData,endIndexData)]);
	  
	  log_prob = std::min<double>(0.0, new_prob - old_prob);
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
	for(int j = 0; j < n_samp_max; ++j){
	  proposedIndivHist = newInfHist(indiv,_);
	  indivHist = newInfHist(indiv,_);
	  // Which year under consideration
	  // Get index for the matrix Z (not index for individual's inf hist)
	  // Note that locs[j] is the index in the individual's inf hist
	  year = sample_years(locs[j]);
	  // Get number of individuals that were alive and/or infected in that year,
	  // less the current individual
	  // Number of infections in this year, less infection status of this individual in this year
	  m = sum(newInfHist(_,year)) - newInfHist(indiv,year);      
	  // Number alive is number alive overall, less this individual
	  n = n_alive(year) - 1;
	  // Work out proposal ratio - prior from alpha, beta and number of other infections
	  ratio = (m + alpha)/(n + alpha + beta);
	  // Propose 1 or 0 based on this ratio
	  rand1 = R::runif(0,1);
	  if(rand1 < ratio){
	    new_entry = 1;
	    proposedIndivHist(year) = 1;
	  } else {
	    new_entry = 0;
	    proposedIndivHist(year) = 0;
	  }
	
	  // If proposing a change, need to check likelihood ratio
	  if(new_entry != newInfHist(indiv,year)){
	    old_prob = likelihood_data_individual(pars, indivHist, circulationTimes, circulationMapIndices,
						  samplingTimes[Range(startIndexSamples, endIndexSamples)], 
						  indicesTitreDataSample[Range(startIndexSamples,endIndexSamples)], 
						  measuredMapIndices[Range(startIndexData,endIndexData)],
						  antigenicMapLong, 
						  antigenicMapShort,  
						  n_strains,
						  data[Range(startIndexData,endIndexData)]);
	    new_prob = likelihood_data_individual(pars, proposedIndivHist,
						  circulationTimes, circulationMapIndices,
						  samplingTimes[Range(startIndexSamples, endIndexSamples)], 
						  indicesTitreDataSample[Range(startIndexSamples,endIndexSamples)], 
						  measuredMapIndices[Range(startIndexData,endIndexData)],
						  antigenicMapLong, 
						  antigenicMapShort,  
						  n_strains,
						  data[Range(startIndexData,endIndexData)]);

	  
	    log_prob = std::min<double>(0.0, new_prob - old_prob);
	    rand1 = R::runif(0,1);
	    if(log(rand1) < log_prob){
	      // Update the entry in the new matrix Z
	      newInfHist(indiv, year) = new_entry;
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
			    IntegerVector sampledIndivs, 
			    IntegerVector ageMask,
			    IntegerVector moveSizes, 
			    IntegerVector nInfs,
			    double alpha, 
			    double beta, 
			    NumericVector randNs) {
  // Copy input matrix
  arma::mat newInfHist = infHist;
  IntegerVector locs; // Locations to be updated
  arma::uvec locs1;
  arma::mat x;
  arma::mat y;
  IntegerVector samps;
  
  int index = 0;
  
  int maxI_indiv;
  int maxI = newInfHist.n_cols;
  int indiv;
  int k;
  int nInf;
  int n;
  int moveMax;
  int move;
  int id1;
  int id2;
  int tmp;
  
  double rand1;
  double ratio;
  
  // For each sampled individual
  for(int i = 0; i < sampledIndivs.size(); ++i){

    // Isolate that individual's infection histories
    indiv = sampledIndivs[i]-1;
    //Rcpp::Rcout << "individual: " << indiv << std::endl;
    nInf = nInfs[indiv];
    x = newInfHist.submat(indiv, ageMask[indiv]-1, indiv, maxI-1);
    samps = seq_len(x.n_cols);
    // With 50% probability, add/remove infections or swap infections
    if(randNs[i] < 1.0/2.0){
      /*
      Rcpp::Rcout << "Indiv: " << indiv << std::endl;
      Rcpp::Rcout << "ageMask: " << ageMask[indiv]-1 << std::endl;
      Rcpp::Rcout << "nInf: " << nInf << std::endl;
      Rcpp::Rcout << "Inf hist: " << x << std::endl;
      Rcpp::Rcout << "Samps: " << samps << std::endl;
      */
      // Sample N random locations
      locs = RcppArmadillo::sample(samps, nInf, FALSE, NumericVector::create());
      //Rcpp::Rcout << "Available locs: " << locs << std::endl;
      locs1 = as<arma::uvec>(locs)-1;
      //Rcpp::Rcout << "Locs chosen: " << locs1 << std::endl;
      //y = newInfHist.row(indiv);
      //Rcpp::Rcout << "y: " << y << std::endl;
      y = x.elem(locs1);
      //Rcpp::Rcout << "locations chosen contents: " << y << std::endl;
      // Count the number of 1s and 0s
      k = accu(x) - accu(y);
      //Rcpp::Rcout << "Infections less sampled: " << k << std::endl;
      n = x.size() - nInf;
      //Rcpp::Rcout << "N less sampled: " << n << std::endl;

      
      // For each sampled location, choose to turn into a 1 or 0 depending
      // on the beta binomial distribution.
      for(int j = 0; j < nInf; ++j){
	//Rcpp::Rcout << "Alpha: " << alpha << "; beta: " << beta << std::endl;
        ratio = (alpha + k)/(alpha + beta + n);
	//Rcpp::Rcout << "ratio: " << ratio << std::endl;
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
      for(int j = 0; j < nInf; j++){
        id1 = floor(R::runif(0,1)*x.size());
        moveMax = moveSizes[indiv];
        move = floor(R::runif(0,1)*2*moveMax) - moveMax;
        id2 = id1 + move;
	while(id2 < 0) id2 += maxI_indiv;
        //if(id2 < 0) id2 = (move + id1) % maxI_indiv;
	if(id2 >= maxI_indiv) id2 = (id2 % maxI_indiv);
	tmp = x[id1];
        x[id1] = x[id2];
        x[id2] = tmp;
      }
    }
    newInfHist.submat(indiv, ageMask[indiv]-1, indiv, maxI-1) = x;
  }

  return(newInfHist);
}

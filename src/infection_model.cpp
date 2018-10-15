#include <Rcpp.h>
#include "wane_function.h"
#include "boosting_functions.h"
using namespace Rcpp;

#define MAX(a,b) ((a) < (b) ? (b) : (a)) // define MAX function for use later
  
void setup_waning_and_masked_cumulative(const NumericVector &theta,
					const IntegerVector &infection_history,
					IntegerVector &cumu_infection_history, 
					IntegerVector &masked_infection_history,
					NumericVector &waning,
					const NumericVector &infection_times,
					const int &max_infections, 
					const double &sampling_time,
					const int &wane_type,
					const double DOB=0,
					const Nullable<List> &additional_arguments=R_NilValue){
  double circulation_time=infection_times[0];
  double time = sampling_time - circulation_time;
  double wane = theta["wane"];
  double val;  
  
  /* =========================================
  ======== SET UP WANING RATE VECTOR =========
  ========================================= */ 
  // If not linear
  if(wane_type == 1){
    // Calculate the waning at the time since infection
    val= wane_function(theta, time, wane);
    waning[0] = MAX(0, 1.0-val);
  } else {
    waning[0] = MAX(0, 1.0-wane*time);// Else if linear
  }

  /* =========================================
     ====== SET UP CUMULATIVE INFECTIONS =====
     ========================================= */ 
  if(circulation_time > sampling_time) masked_infection_history[0]=0;
  else masked_infection_history[0] = infection_history[0];

  cumu_infection_history[0] = masked_infection_history[0];

  for(int i = 1; i < max_infections; ++i){
    circulation_time = infection_times[i];
    time = (sampling_time - circulation_time);

    /* =========================================
       ====== SET UP CUMULATIVE INFECTIONS =====
       ========================================= */
    /* Check if isolation time is after the sampling time.
       if so, then we do not test against this strain */ 
    if(circulation_time > sampling_time){
      masked_infection_history[i] = 0;
    } else {
      masked_infection_history[i] = infection_history[i];
    }
    cumu_infection_history[i] = masked_infection_history[i] + cumu_infection_history[i-1];

    /* =========================================
       ======== SET UP WANING RATE VECTOR =========
       ========================================= */ 
    // If not linear
    if(wane_type == 1){
      val= wane_function(theta, time, wane);
      waning[i] = MAX(0, 1.0-val);
    } else {
      waning[i] = MAX(0, 1.0-wane*time); // Else linear
    }
  }  
}




//' Model function
//'
//' The main model solving function for a single individual.
//' NOTES:
//' - Do we want infection history to be a vector of infection times?
//' - Treat the contents of infectionHistory as a parameter (ie. exposure type)
//' @param theta NumericVector, the named vector of model parameters
//' @param infectionHistory IntegerVector, the vector of 1s and 0s showing presence/absence of infection for each possible time. 
//' @param infectionTimes NumericVector, the actual times of circulation that the infection history vector corresponds to
//' @param infectionMapIndices IntegerVector, which entry in the melted antigenic map that these infection times correspond to
//' @param samplingTime double, the real time that the sample was taken
//' @param measurementMapIndices IntegerVector, the indices of all measured strains in the melted antigenic map
//' @param antigenicMapLong NumericVector, the collapsed cross reactivity map for long term boosting, after multiplying by sigma1
//' @param antigenicMapShort NumericVector, the collapsed cross reactivity map for short term boosting, after multiplying by sigma2
//' @param numberStrains int, the maximum number of infections that an individual could experience
//' @return NumericVector of predicted titres for each entry in measurementMapIndices
//' @useDynLib serosolver
//' @export
//[[Rcpp::export]]
NumericVector infection_model_indiv(const NumericVector &theta, // Parameter vector
				    const IntegerVector &infectionHistory, // vector of 1s and 0s for infections
				    const NumericVector &infectionTimes, // Time of these infections
				    const IntegerVector &infectionMapIndices, // Where these infection times fit in antigenic map
				    const double &samplingTime,  // This sampling time
				    const IntegerVector &measurementMapIndices, // Indices of measured strains in antigenic map
				    const NumericVector &antigenicMapLong,
				    const NumericVector &antigenicMapShort, 
				    const int &numberStrains, // Maximum number of infections an individual could experience, if alive the whole time
				    const double &DOB=0,
				    const Nullable<List> &additional_arguments=R_NilValue
				    ){
  double circulation_time;
  
  // Which waning function type should we use
  // 0 is linear decrease
  // 1 is piecewise linear
  int waneType = theta["wane_type"]; 
  // If titre dependent boosting or not
  bool titre_dependent_boosting = theta["titre_dependent"] == 1;

  // We will need to loop over each strain that was tested
  int n_samples = measurementMapIndices.size(); // Number of time points sampled
  int max_infections = infectionTimes.size(); // max number of infections is one for each strain
  
  // Only recording titres for which we have data
  NumericVector predictedTitre(n_samples);
  NumericVector monitored_titres(max_infections);
  // But need to record infection info for each strain that individual could have been
  // infected with
  IntegerVector cumInfectionHistory(max_infections);
  IntegerVector maskedInfectionHistory(max_infections); // To avoid solving the model for infections that didn't happen
  NumericVector waning(max_infections);

  // Set up cumulative infection history, masked infection history and waning vector
  setup_waning_and_masked_cumulative(theta,
				     infectionHistory,
				     cumInfectionHistory, 
				     maskedInfectionHistory, 
				     waning,				     
				     infectionTimes, 
				     max_infections, 
				     samplingTime,
				     waneType, 
				     DOB,
				     additional_arguments);
  add_multiple_infections_boost(predictedTitre, monitored_titres,
				theta, 
				infectionTimes, cumInfectionHistory, maskedInfectionHistory, 
				infectionMapIndices, measurementMapIndices, 
				antigenicMapLong, antigenicMapShort,
				waning, 
				numberStrains, 
				titre_dependent_boosting, 
				DOB,
				additional_arguments
				);
  return(predictedTitre);
}


//' @export
//[[Rcpp::export]]
NumericVector titre_data_individual(const NumericVector &theta, 
				    const IntegerVector &infectionHistory, 
				    const NumericVector &circulationTimes, 
				    const IntegerVector &circulationMapIndices,
				    const NumericVector &samplingTimes,
				    const IntegerVector &dataIndices,
				    const IntegerVector &measuredMapIndices, 
				    const NumericVector &antigenicMapLong, 
				    const NumericVector &antigenicMapShort,
				    const int &numberStrains,
				    const double &DOB=0,
				    const Nullable<List>& additional_arguments=R_NilValue){
  int numberSamples = samplingTimes.size();
  int numberMeasuredStrains = measuredMapIndices.size();
  NumericVector titres(numberMeasuredStrains);

  int startIndex = 0;
  int endIndex = 0;

  LogicalVector indices = infectionHistory > 0;

  IntegerVector conciseInfHist = infectionHistory[indices];
  NumericVector infectionTimes = circulationTimes[indices];
  IntegerVector infMapIndices = circulationMapIndices[indices];

  for(int i = 0; i < numberSamples; ++i){
    endIndex = startIndex + dataIndices[i] - 1;
    titres[Range(startIndex, endIndex)] = infection_model_indiv(theta,conciseInfHist,infectionTimes,infMapIndices,
								samplingTimes[i],measuredMapIndices[Range(startIndex,endIndex)],
								antigenicMapLong, antigenicMapShort,numberStrains,
								DOB, additional_arguments);
    startIndex = endIndex + 1;
  }
  return(titres);
}

//' @export
//[[Rcpp::export]]
NumericVector titre_data_group(const NumericVector &theta, 
			       const IntegerMatrix &infectionHistories, 
			       const NumericVector &circulationTimes,
			       const IntegerVector &circulationMapIndices,
			       const NumericVector &samplingTimes,
			       const IntegerVector &indicesTitreDataSample, // How many rows in titre data correspond to each individual, sample and repeat?
			       const IntegerVector &indicesTitreDataOverall, // How many rows in the titre data correspond to each individual?
			       const IntegerVector &indicesSamples, // Split the sample times and runs for each individual
			       const IntegerVector &measuredMapIndices, // For each titre measurement, corresponding entry in antigenic map
			       const NumericVector &antigenicMapLong, 
			       const NumericVector &antigenicMapShort,
			       const NumericVector &DOBs,
			       const Nullable<List> &additional_arguments=R_NilValue
			       ){
  int n = infectionHistories.nrow();
  int n_strains = infectionHistories.ncol();

  NumericVector titres(measuredMapIndices.size());
  
  bool check_additional_arguments = additional_arguments.isNotNull(); // Precompute this check so only have to do it once
  bool titre_dependent_boosting = theta["titre_dependent"] == 1;

  // These indices allow us to step through the titre data vector
  // as if it were a matrix ie. number of rows for each individual
  // at a time
  int startIndexSamples;
  int endIndexSamples;
  int startIndexData;
  int endIndexData;
  int DOB;

  for(int i=1; i <= n; ++i){
    startIndexSamples = indicesSamples[i-1];
    endIndexSamples = indicesSamples[i] - 1;

    startIndexData = indicesTitreDataOverall[i-1];
    endIndexData = indicesTitreDataOverall[i] - 1;

    DOB = DOBs[i-1];

    titres[Range(startIndexData, endIndexData)] = titre_data_individual(theta,   // Vector of named model parameters
									infectionHistories(i-1,_), // Vector of infection history for individual i
									circulationTimes, // Vector of all virus circulation times, same length and ncol infectionHistories
									circulationMapIndices, // Gives the corresponding index in the antigenic map vector
									samplingTimes[Range(startIndexSamples, endIndexSamples)],  // Get sampling times for this individual
									indicesTitreDataSample[Range(startIndexSamples,endIndexSamples)], // The 
									measuredMapIndices[Range(startIndexData,endIndexData)],  // For the indices in the antigenic map to which each titre corresponds
									antigenicMapLong, 
									antigenicMapShort,  
									n_strains,
									DOB, 
									additional_arguments); // The total number of strains that circulated
  }

  return(titres);
}



//' Calculate likelihood
//'
//' Calculates the likelihood of a given set of observed titres given predicted titres. Based on truncated discritised normal.
//' @param expected NumericVector, as returned by \code{\link{infection_model_indiv}}
//' @param data NumericVector, the vector of observed titres
//' @param theta NumericVector, the vector of named model parameters, requiring MAX_TITRE and error
//' @param titre_shifts NumericVector, OPTIONAL if using measurement bias, gives the shift to add to each expected titre
//' @return a single log likelihood
//' @export
//[[Rcpp::export]]
double likelihood_titre(const NumericVector &expected, 
			const NumericVector &data, 
			const NumericVector &theta,
			const Nullable<NumericVector> &titre_shifts=R_NilValue){
  NumericVector expected1 = expected;
  if(titre_shifts.isNotNull()){
    NumericVector _titre_shifts(titre_shifts);
    expected1 = expected1 + _titre_shifts;
  }
  int n = expected.size();
  double lnlike = 0;
  
  for(int i=0; i < n; ++i){
    if(!NumericVector::is_na(data[i])){
      if(data[i] > theta["MAX_TITRE"]){
	lnlike += R::pnorm(theta["MAX_TITRE"], expected[i], theta["error"],0,1);
      } else if(data[i] <= 0){
	lnlike += R::pnorm(1, expected[i], theta["error"],1,1);
      } else {
	lnlike += log(R::pnorm(data[i]+1, expected[i],theta["error"], 1,0) - 
		      R::pnorm(data[i], expected[i], theta["error"],1,0));
      }
    }
  }
  return(lnlike);
}



//' @export
//[[Rcpp::export]]
double likelihood_data_individual(const NumericVector &theta, 
				  const IntegerVector &infectionHistory, 
				  const NumericVector &circulationTimes, 
				  const IntegerVector &circulationMapIndices,
				  const NumericVector &samplingTimes,
				  const IntegerVector &dataIndices,
				  const IntegerVector &measuredMapIndices, 
				  const NumericVector &antigenicMapLong, 
				  const NumericVector &antigenicMapShort,
				  const int &numberStrains,
				  const NumericVector &data,
				  const Nullable<NumericVector> &titre_shifts,
				  const double &DOB=0,
				  const Nullable<List> &additional_arguments=R_NilValue
				  ){
  int numberSamples = samplingTimes.size();
  int numberMeasuredStrains = measuredMapIndices.size();
  double lnlike=0;
  NumericVector titres(numberMeasuredStrains);

  int startIndex = 0;
  int endIndex = 0;

  LogicalVector indices = infectionHistory > 0;

  IntegerVector conciseInfHist = infectionHistory[indices];
  NumericVector infectionTimes = circulationTimes[indices];
  IntegerVector infMapIndices = circulationMapIndices[indices];

  for(int i = 0; i < numberSamples; ++i){ 
    endIndex = startIndex + dataIndices[i] - 1;
    titres[Range(startIndex, endIndex)] = infection_model_indiv(theta,
								conciseInfHist,
								infectionTimes,
								infMapIndices,
								samplingTimes[i],
								measuredMapIndices[Range(startIndex,endIndex)],
								antigenicMapLong, antigenicMapShort,
								numberStrains, DOB, additional_arguments);
    startIndex = endIndex + 1;
  }
  lnlike = likelihood_titre(titres, data, theta, titre_shifts);
  return(lnlike);
}


//' Sums a vector based on bucket sizes
//'
//' Given a vector (a) and another vector of bucket sizes, returns the summed vector (a)
//' @param a the vector to be bucketed
//' @param buckets the vector of bucket sizes to sum a over
//' @return the vector of summed a
//' @export
//[[Rcpp::export]]
NumericVector sum_buckets(NumericVector a, NumericVector buckets){
  NumericVector results(buckets.size());
  int index = 0;
  for(int i = 0; i < buckets.size(); ++i){
    results[i] = 0;
    for(int j = 0; (j < buckets[i]) & (index < a.size()); ++j){
      results[i] += a[index++];
    }
  }
  return(results);
}

//' Convert melted antigenic map to cross reactivity
//'
//' Multiplies all elements of the provided vector, x such that y = 1 - sigma*x. Also makes sure that no calculated value is less than 0
//' @param x the melted antigenic map
//' @param sigma the cross reactivity waning parameter
//' @return a vector of cross reactivity
// [[Rcpp::export]]
NumericVector create_cross_reactivity_vector(NumericVector x, double sigma) {
  NumericVector x2(x.size());
  for(int i = 0; i < x2.size(); ++i){
    x2[i] = MAX(1 - x[i]*sigma,0);
  }
  return(x2);
}


 
//' Original model reimplementation
//' @export
// [[Rcpp::export]]
NumericVector c_model_original(int n, // max number of infections
			       int nsamp, // number of samples tested against (ie. number of titres)
			       IntegerVector x, // infection history vector
			       NumericVector theta,
			       NumericVector dd, // Long term CR matrix
			       NumericVector dd2, // short term CR matrix
			       int t_sample){ // index of time that sample was taken (ie. 1 = 1968)
  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /* Calculate lambda */	
  double T_2 = theta["tau"];
  double wane = theta["wane"];
  double mu = theta["mu"];
  double mu2 = theta["mu_short"]; 
  // as sigma = theta[4]
	
  NumericVector titrepred(nsamp);
  IntegerVector maskedInfectionHistory(n);
  NumericVector distanceFromTest(n);
  IntegerVector cumInfectionHistory(n);
  NumericVector x1(n);
	
#define MAX(a,b) ((a) < (b) ? (b) : (a)) // define MAX function for use later
  
  /* Add for loop over k*/
	
  int k;
  int i;
  int j;
  int m;
	
  double xx2; 
  // For each strain tested
  for (k=0; k<nsamp; k++){ // Iterate over samples tested against

    xx2=0;

    // Make a masked infection history
    // Was the individual infected in each year?
    j=(t_sample-1); // fix test year. Need (t-1) as index from 0
    // Iterate through all possible infection years
    for (m=0;m<n;m++) {
      if (m <= j) {
	maskedInfectionHistory[m] = x[m];
      } else {
	maskedInfectionHistory[m] = 0;
      }
    }
    
    // Make an index for waning
    // For all possible infection years
    // how far is the tested year from a year of infection?
    for (m=0;m<n;m++) {
      // distanceFromTest[m] = exp(-wane * (j-m )); // Distance from test year 
      distanceFromTest[m] = MAX(0, 1 - wane * (j-m) ); // Distance from test year
    }
    
    // Make a cumulative infection history
    cumInfectionHistory[0] = maskedInfectionHistory[0];
    for (m=1;m<n;m++) {
      cumInfectionHistory[m] = cumInfectionHistory[m-1] + 
	maskedInfectionHistory[m];
    }
    
    /* Calculate expected titre	- note k indexed from 0 */
    /* Note that waning is currently linked with back boosting */
    /* dd is long term cross-reaction, dd2 is short-term */
    
    for (i=0; i<n; i++){
      x1[i] = maskedInfectionHistory[i] *
	// exp(-1.0 * T_2 * ( cumInfectionHistory[i]  - 1.0)) *
	MAX(0, 1.0 - T_2 * ( cumInfectionHistory[i]  - 1.0)) * // Antigenic seniority
	//(pow(1.0 + T_1 , (total_inf - cumInfectionHistory[i])) ) * REMOVED Tau 1
	(mu * dd[k*n+i] + mu2 * dd2[k*n+i] * distanceFromTest[i] );
    }
    // Sum up contribution of all infections
    for (i=0; i<n; i++){
      xx2 =  xx2 + x1[i];
    }
    
    titrepred[k] = xx2;
    
    
  } // end sample loop (k)
  return(titrepred);
}

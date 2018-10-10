#include <Rcpp.h>
using namespace Rcpp;

#define MAX(a,b) ((a) < (b) ? (b) : (a)) // define MAX function for use later
  
double likelihood_titre(NumericVector expected, NumericVector data, NumericVector theta, Nullable<NumericVector> to_add);

//' Model function
//'
//' The main model solving function for a single individual.
//' NOTES:
//' - Do we want infection history to be a vector of infection times?
//' - Treat the contents of infectionHistory as a parameter (ie. exposure type)
//' @param theta NumericVector, the named vector of model parameters
//' @param mus NumericVector, the vector of boosting parameters corresponding to each circulation time
//' @param boostingVecIndices, which index in the mus vector does each potential infection correspond to?
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
NumericVector infection_model_indiv_mus(NumericVector theta, // Parameter vector
					NumericVector mus,
					IntegerVector infectionHistory, // vector of 1s and 0s for infections
					NumericVector infectionTimes, // Time of these infections
					IntegerVector boostingVecIndices,
					IntegerVector infectionMapIndices, // Where these infection times fit in antigenic map
					double samplingTime,  // This sampling time
					IntegerVector measurementMapIndices, // Indices of measured strains in antigenic map
					NumericVector antigenicMapLong,
					NumericVector antigenicMapShort, 
					int numberStrains // Maximum number of infections an individual could experience
					){
  // Extract model parameters
  double mu = 0;
  double mu_short = theta["mu_short"];
  double tau = theta["tau"];
  double wane = theta["wane"];
  
  // We will need to loop over each strain that was tested
  int n_samples = measurementMapIndices.size(); // Number of time points sampled
  int max_infections = infectionTimes.size(); // max number of infections is one for each strain
  double tmpTitre=0;
  
  // Only recording titres for which we have data
  NumericVector predictedTitre(n_samples);
  
  // But need to record infection info for each strain that individual could have been
  // infected with
  IntegerVector cumInfectionHistory(max_infections);
  IntegerVector maskedInfectionHistory(max_infections);
  NumericVector waning(max_infections);

  double circulation_time;

  /*
  Rcpp::Rcout << "Number of titres performed: " << n_samples << std::endl;
  Rcpp::Rcout << "Maximum number of infections: " << max_infections << std::endl;
  Rcpp::Rcout << "Number of strains: " << numberStrains << std::endl;
  Rcpp::Rcout << "Waning rate: " << wane << std::endl;
  */
  // Set up cumulative infection history and waning vector
  /* Check if isolation time is after the sampling time.
     if so, then we do not test against this strain */ 
  circulation_time = infectionTimes[0];
  if(circulation_time > samplingTime) maskedInfectionHistory[0]=0;
  else maskedInfectionHistory[0] = infectionHistory[0];
  waning[0] = MAX(0, 1.0-wane*(samplingTime - circulation_time));
  cumInfectionHistory[0] = maskedInfectionHistory[0];

  /* For strains that circulated before the isolation time,
     add to cumulative infection history */
  for(int i=1; i < max_infections; ++i){
    /* At which time did this strain circulate?
       This might change to a range of times in the future */
    circulation_time = infectionTimes[i];
    // If circulated after the sample, then could not have been infected
    if(circulation_time > samplingTime){
      maskedInfectionHistory[i]=0;
    } else {
      maskedInfectionHistory[i] = infectionHistory[i];
    }

    // Make cumulative infection history
    cumInfectionHistory[i] = maskedInfectionHistory[i] + cumInfectionHistory[i-1];   

    /* Get waning rate for this strain
       Note that circulation_time is implicitly the time at which
       an individual would have been infected */
    waning[i] = MAX(0, 1.0-wane*(samplingTime-circulation_time));
  }
  /*
  Rcpp::Rcout << "Cumulative infecion history: " << cumInfectionHistory << std::endl;
  Rcpp::Rcout << "Infection map indices: " <<  infectionMapIndices << std::endl;
  Rcpp::Rcout << "Infection times: " <<  infectionTimes << std::endl;
  Rcpp::Rcout << "Infection history: " <<  infectionHistory << std::endl;
  Rcpp::Rcout << "Measurement map indices: " <<  measurementMapIndices << std::endl;
  Rcpp::Rcout << "Mus: " <<  mus << std::endl;
  Rcpp::Rcout << "Mu indices: " <<  boostingVecIndices << std::endl;
  */


  // For each strain we are testing against, find predicted titre
  for(int k=0; k < n_samples; ++k){
    tmpTitre=0;
    // Sum contributions from all infections
    // Note that measurementMapIndices[k] finds the correct entry in the
    // antigenicMap vector for the *tested* strain (row), whereas
    // infectionMapIndices[i] finds the entry for the *infecting* strain (column)
    for(int i=0; i < max_infections; ++i){
      mu = mus[boostingVecIndices[infectionMapIndices[i]]];
      /*
	Rcpp::Rcout << "Sample no: " << k << std::endl;
	Rcpp::Rcout << "Infection no: " << i << std::endl;
	Rcpp::Rcout << "Infection map index: " << infectionMapIndices[i] << std::endl;
	Rcpp::Rcout << "Boosting vector index: " << boostingVecIndices[infectionMapIndices[i]] << std::endl;
	Rcpp::Rcout << "Realised mu: " << mu << std::endl;
      */
      ///////////////////////////////
      // THE ACTUAL MODEL
      tmpTitre += maskedInfectionHistory[i]* // Ignore infections that couldn't have happened
	MAX(0, 1.0 - tau*(cumInfectionHistory[i] - 1.0))* // Antigenic seniority
	(mu*antigenicMapLong[measurementMapIndices[k]*numberStrains+infectionMapIndices[i]] +  // Long term boost
	 mu_short*antigenicMapShort[measurementMapIndices[k]*numberStrains+infectionMapIndices[i]]* // Short term cross reactive boost
	 waning[i]); // Waning rate
      ////////////////////////////
    }
    predictedTitre[k] = tmpTitre;
  }
  return(predictedTitre);
}


//' @export
//[[Rcpp::export]]
NumericVector titre_data_individual_mus(NumericVector theta, 
					NumericVector mus,
					IntegerVector infectionHistory, 
					NumericVector circulationTimes, 
					IntegerVector circulationMapIndices,
					IntegerVector musIndices,
					NumericVector samplingTimes,
					IntegerVector dataIndices,
					IntegerVector measuredMapIndices, 
					NumericVector antigenicMapLong, 
					NumericVector antigenicMapShort,
					int numberStrains
					){
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
    //Rcpp::Rcout << "Mu indices: " << musIndices << std::endl;
    //Rcpp::Rcout << "Mus: " << mus << std::endl;
    //Rcpp::Rcout << "Infection times: " << infectionTimes << std::endl;
    titres[Range(startIndex, endIndex)] = infection_model_indiv_mus(theta,mus, conciseInfHist,infectionTimes,
								musIndices, infMapIndices,samplingTimes[i],
								measuredMapIndices[Range(startIndex,endIndex)],
								antigenicMapLong, antigenicMapShort,numberStrains);
    startIndex = endIndex + 1;
  }
  return(titres);
}

// @export
//[[Rcpp::export]]
NumericVector titre_data_group_mus(NumericVector theta, 
				   NumericVector mus,
				   IntegerMatrix infectionHistories, 
				   NumericVector circulationTimes,
				   IntegerVector circulationMapIndices,
				   IntegerVector musIndices,
				   NumericVector samplingTimes,
				   IntegerVector indicesTitreDataSample, // How many rows in titre data correspond to each individual, sample and repeat?
				   IntegerVector indicesTitreDataOverall, // How many rows in the titre data correspond to each individual?
				   IntegerVector indicesSamples, // Split the sample times and runs for each individual
				   IntegerVector measuredMapIndices, // For each titre measurement, corresponding entry in antigenic map
				   NumericVector antigenicMapLong, 
				   NumericVector antigenicMapShort){
  int n = infectionHistories.nrow();
  int n_strains = infectionHistories.ncol();

  NumericVector titres(measuredMapIndices.size());
  
  // These indices allow us to step through the titre data vector
  // as if it were a matrix ie. number of rows for each individual
  // at a time
  int startIndexSamples;
  int endIndexSamples;
  int startIndexData;
  int endIndexData;

  for(int i=1; i <= n; ++i){
    //Rcpp::Rcout << "Individual: " << i << std::endl;
    startIndexSamples = indicesSamples[i-1];
    endIndexSamples = indicesSamples[i] - 1;

    startIndexData = indicesTitreDataOverall[i-1];
    endIndexData = indicesTitreDataOverall[i] - 1;
    titres[Range(startIndexData, endIndexData)] = titre_data_individual_mus(theta,   // Vector of named model parameters
									mus,
									infectionHistories(i-1,_), // Vector of infection history for individual i
									circulationTimes, // Vector of all virus circulation times, same length and ncol infectionHistories
									circulationMapIndices, // Gives the corresponding index in the antigenic map vector
									musIndices,
									samplingTimes[Range(startIndexSamples, endIndexSamples)],  // Get sampling times for this individual
									indicesTitreDataSample[Range(startIndexSamples,endIndexSamples)], // The 
									measuredMapIndices[Range(startIndexData,endIndexData)],  // For the indices in the antigenic map to which each titre corresponds
									antigenicMapLong, 
									antigenicMapShort,  
									n_strains); // The total number of strains that circulated
  }

  return(titres);
}


//' @export
//[[Rcpp::export]]
double likelihood_data_individual_mus(NumericVector theta, 
				      NumericVector mus,
				      IntegerVector infectionHistory, 
				      NumericVector circulationTimes, 
				      IntegerVector circulationMapIndices,
				      IntegerVector musIndices,
				      NumericVector samplingTimes,
				      IntegerVector dataIndices,
				      IntegerVector measuredMapIndices, 
				      NumericVector antigenicMapLong, 
				      NumericVector antigenicMapShort,
				      int numberStrains,
				      NumericVector data
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
    titres[Range(startIndex, endIndex)] = infection_model_indiv_mus(theta,mus,
								    conciseInfHist,
								    infectionTimes,
								    musIndices,
								    infMapIndices,
								    samplingTimes[i],
								    measuredMapIndices[Range(startIndex,endIndex)],
								    antigenicMapLong, antigenicMapShort,
								    numberStrains);
    startIndex = endIndex + 1;
  }
  lnlike = likelihood_titre(titres, data, theta, R_NilValue);
  return(lnlike);
}

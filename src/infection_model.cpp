#include <Rcpp.h>
using namespace Rcpp;

#define MAX(a,b) ((a) < (b) ? (b) : (a)) // define MAX function for use later
  
//' Model function
//'
//' The main model solving function for a single individual.
//' NOTES:
//' - Do we want infection history to be a vector of infection times?
//' - Treat the contents of infectionHistory as a parameter (ie. exposure type)
//' @param theta NumericVector, the named vector of model parameters
//' @param infectionHistory NumericVector, the vector of 1s and 0s showing presence/absence of infection for each possible time. 
//' @param samplingTime double, the real time that the sample was taken
//' @param strainIsolationTimes NumericVector, the vector of times at which each virus strain circulated
//' @param antigenicMapLong NumericVector, the collapsed cross reactivity map for long term boosting, after multiplying by sigma1
//' @param antigenicMapShort NumericVector, the collapsed cross reactivity map for short term boosting, after multiplying by sigma2
//' @return NumericVector of predicted titres for each strainIsolationTime
//' @useDynLib serosolver
//' @export
//[[Rcpp::export]]
NumericVector infection_model_indiv(NumericVector theta, NumericVector infectionHistory,
				    double samplingTime, NumericVector strainIsolationTimes, IntegerVector virusIndices,
				    NumericVector antigenicMapLong,NumericVector antigenicMapShort, int numberStrains){
  // Extract model parameters
  double mu = theta["mu"];
  double mu_short = theta["mu_short"];
  double tau = theta["tau"];
  double wane = theta["wane"];

  // We will need to loop over each strain that was tested
  int n_samples = strainIsolationTimes.size(); // Number of time points sampled
  int max_infections = numberStrains; // max number of infections is one for each strain
  double tmpTitre=0;

  
// Only recording titres for which we have data
  NumericVector predictedTitre(n_samples);
  
  // But need to record infection info for each strain that individual could have been
  // infected with
  NumericVector cumInfectionHistory(max_infections);
  NumericVector maskedInfectionHistory(max_infections);
  NumericVector waning(max_infections);

  double circulation_time;

  // Set up cumulative infection history and waning vector
  /* Check if isolation time is after the sampling time.
     if so, then we do not test against this strain */ 
  circulation_time = strainIsolationTimes[0];
  if(circulation_time > samplingTime) maskedInfectionHistory[0]=0;
  else maskedInfectionHistory[0] = infectionHistory[0];
  waning[0] = MAX(0, 1.0-wane*(samplingTime - circulation_time));
  cumInfectionHistory[0] = maskedInfectionHistory[0];

  /* For strains that circulated before the isolation time,
     add to cumulative infection history */
  //Rcpp::Rcout << "Sampling time: " << samplingTime << std::endl;
  for(int i=1; i < numberStrains; ++i){
    /* At which time did this strain circulate?
       This might change to a range of times in the future */
    circulation_time = strainIsolationTimes[i];
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
  
  // For each strain we are testing against, find predicted titre
  for(int k=0; k < n_samples; ++k){
    tmpTitre=0;
    // Sum contributions from all infections
    // Note that virusIndices[k] finds the correct entry in the
    // antigenicMap vector for the *tested* strain, whereas
    // i finds the entry for the *infecting* strain
    for(int i=0; i < max_infections; ++i){
      ///////////////////////////////
      // THE ACTUAL MODEL
      tmpTitre += maskedInfectionHistory[i]* // Ignore infections that couldn't have happened
	            MAX(0, 1.0 - tau*(cumInfectionHistory[i] - 1.0))* // Antigenic seniority
	            (mu*antigenicMapLong[virusIndices[k]*numberStrains+i] +  // Long term boost
	            mu_short*antigenicMapShort[virusIndices[k]*numberStrains+i]* // Short term cross reactive boost
	            waning[i]); // Waning rate
      ////////////////////////////
    }
    predictedTitre[k] = tmpTitre;
  }
  return(predictedTitre);
}

//' Calculate likelihood
//'
//' Calculates the likelihood of a given set of observed titres given predicted titres. Based on truncated discritised normal.
//' @param expected NumericVector, as returned by \code{\link{infection_model_indiv}}
//' @param data NumericVector, the vector of observed titres
//' @param theta NumericVector, the vector of named model parameters, requiring MAX_TITRE and error
//' @return a single log likelihood
//' @export
//[[Rcpp::export]]
double likelihood_titre(NumericVector expected, NumericVector data, NumericVector theta){
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

//' Individual likelihood
//'
//' Uses \code{\link{infection_model_indiv}} for a single individual and then calls \code{\link{likelihood_titre}}. Works with NA
//' @param theta NumericVector, the named vector of model parameters
//' @param infectionHistory NumericVector, the vector of 1s and 0s showing presence/absence of infection for each possible time. 
//' @param samplingTimes NumericVector, the vector of real times that samples were taken
//' @param strainIsolationTimes NumericVector, the vector of times at which each virus strain circulated
//' @param antigenicMapLong NumericVector, the collapsed cross reactivity map for long term boosting, after multiplying by sigma1
//' @param antigenicMapShort NumericVector, the collapsed cross reactivity map for short term boosting, after multiplying by sigma2
//' @param titres NumericVector, the vector of observed titres.
//' @return a single log likelihood
//' @useDynLib serosolver
//' @export
//[[Rcpp::export]]
double individual_likelihood(NumericVector theta, NumericVector infectionHistory,
			     NumericVector samplingTimes, IntegerVector indivIndices, 
			     NumericVector strainIsolationTimes, IntegerVector indivVirusIndices,
			     NumericVector antigenicMapLong, NumericVector antigenicMapShort, 
			     NumericVector titres, int numberStrains){
  
  int numberSamples = samplingTimes.size();
  NumericVector predictedTitres;
  NumericVector tmpTitres;
  double lnlike = 0;

  // These indices allow us to step through a vector as if it were a 
  // matrix (ie. numberStrains at a time)
  int startIndex = 0;
  int endIndex = indivIndices[0] - 1;
  
  for(int i=0; i < numberSamples; ++i){
    predictedTitres = infection_model_indiv(theta, infectionHistory, samplingTimes[i],
					                                  strainIsolationTimes[Range(startIndex, endIndex)], 
                                            indivVirusIndices[Range(startIndex, endIndex)],
                                            antigenicMapLong, antigenicMapShort, numberStrains);
    tmpTitres = titres[Range(startIndex,endIndex)];
    lnlike += likelihood_titre(predictedTitres, tmpTitres,theta);

    // Increase indices by max number of strains
    startIndex = endIndex+1;
    endIndex += indivIndices[i];
  }
  return(lnlike);
}

//' Group likelihood
//'
//' Uses \code{\link{individual_likelihood}} for each individual and returns a vector of log likelihoods for each individual
//' @param theta NumericVector, the named vector of model parameters
//' @param infectionHistories NumericMatrix, the matrix of 1s and 0s showing presence/absence of infection for each possible time for each individual 
//' @param indicesSampling IntegerVector, the range of indices of the titre vector that correspond to each individual (eg. first 5 titres are for indiv 1, indicesA = c(0,5,...)
//' @param indicesData IntegerVector, the range of indices of the sample vector that correspond to each individual (eg. first 2 sampling times are for indiv 1, indicesB = c(0,2,...)
//' @param samplingTimes NumericVector, the vector of real times that samples were taken
//' @param strainIsolationTimes NumericVector, the vector of times at which each virus strain circulated
//' @param antigenicMapLong NumericVector, the collapsed cross reactivity map for long term boosting, after multiplying by sigma1
//' @param antigenicMapShort NumericVector, the collapsed cross reactivity map for short term boosting, after multiplying by sigma2
//' @param titres NumericVector, the vector of observed titres for all individuals.
//' @param n_strains int, the maximum number of strains that could be tested against
//' @return a NumericVector of log likelihoods for each individual
//' @useDynLib serosolver
//' @export
//[[Rcpp::export]]
NumericVector group_likelihood_vector(NumericVector theta, NumericMatrix infectionHistories, 
				      IntegerVector indicesSamples, IntegerVector indicesData, IntegerVector indicesDataOverall,
				      NumericVector samplingTimes, NumericVector strainIsolationTimes, IntegerVector virusIndices,
				      NumericVector antigenicMapLong, NumericVector antigenicMapShort, 
				      NumericVector titres){
  int n = infectionHistories.nrow();
  int n_strains = infectionHistories.ncol();
  NumericVector lnlikes(n);
  
  // These indices allow us to step through the titre data vector
  // as if it were a matrix ie. number of rows for each individual
  // at a time
  int startIndexSamples;
  int endIndexSamples;
  int startIndexData;
  int endIndexData;

  for(int i=0; i < n; ++i){
    startIndexSamples = indicesSamples[i];
    endIndexSamples = indicesSamples[i+1];

    startIndexData = indicesDataOverall[i];
    endIndexData = indicesDataOverall[i+1] - 1;
    
    endIndexSamples -= 1;
    
    lnlikes[i] = individual_likelihood(theta, 
                                      infectionHistories(i,_),
                                      samplingTimes[Range(startIndexSamples, endIndexSamples)],
                                                   indicesData[Range(startIndexSamples,endIndexSamples)],
                                    
                                      strainIsolationTimes[Range(startIndexData,endIndexData)], 
                                                          virusIndices[Range(startIndexData,endIndexData)],
                                      antigenicMapLong, 
                                      antigenicMapShort, 
                                      titres[Range(startIndexData,endIndexData)],
                                            n_strains);
  }
  return(lnlikes);  
}


//' Group likelihood
//'
//' Uses \code{\link{individual_likelihood}} for each individual and returns a summed log likelihood for all individuals
//' @param theta NumericVector, the named vector of model parameters
//' @param infectionHistories NumericMatrix, the matrix of 1s and 0s showing presence/absence of infection for each possible time for each individual 
//' @param indicesA IntegerVector, the range of indices of the titre vector that correspond to each individual (eg. first 5 titres are for indiv 1, indicesA = c(0,5,...)
//' @param indicesB IntegerVector, the range of indices of the sample vector that correspond to each individual (eg. first 2 sampling times are for indiv 1, indicesB = c(0,2,...)
//' @param samplingTimes NumericVector, the vector of real times that samples were taken
//' @param strainIsolationTimes NumericVector, the vector of times at which each virus strain circulated
//' @param antigenicMapLong NumericVector, the collapsed cross reactivity map for long term boosting, after multiplying by sigma1
//' @param antigenicMapShort NumericVector, the collapsed cross reactivity map for short term boosting, after multiplying by sigma2
//' @param titres NumericVector, the vector of observed titres for all individuals.
//' @return a single log likelihood
//' @useDynLib serosolver
//' @export
//[[Rcpp::export]]
double group_likelihood_total(NumericVector theta, NumericMatrix infectionHistories, 
                                      IntegerVector indicesSamples, IntegerVector indicesData, IntegerVector indicesDataOverall,
                                      NumericVector samplingTimes, NumericVector strainIsolationTimes, IntegerVector virusIndices,
                                      NumericVector antigenicMapLong, NumericVector antigenicMapShort, 
                                      NumericVector titres){
  int n = infectionHistories.nrow();
  int n_strains = infectionHistories.ncol();
  double lnlike=0;
  
  // These indices allow us to step through the titre data vector
  // as if it were a matrix ie. number of rows for each individual
  // at a time
  int startIndexSamples;
  int endIndexSamples;
  int startIndexData;
  int endIndexData;
  
  for(int i=0; i < n; ++i){
    startIndexSamples = indicesSamples[i];
    endIndexSamples = indicesSamples[i+1];
    
    startIndexData = indicesDataOverall[i];
    endIndexData = indicesDataOverall[i+1] - 1;
    
    endIndexSamples -= 1;
    
    lnlike += individual_likelihood(theta, 
                                       infectionHistories(i,_),
                                       samplingTimes[Range(startIndexSamples, endIndexSamples)],
                                                    indicesData[Range(startIndexSamples,endIndexSamples)],
                                                               strainIsolationTimes[Range(startIndexData,endIndexData)], 
                                                                                   virusIndices[Range(startIndexData,endIndexData)],
                                                                                               antigenicMapLong, 
                                                                                               antigenicMapShort, 
                                                                                               titres[Range(startIndexData,endIndexData)],
                                                                                                     n_strains);
  }
  return(lnlike);  
}

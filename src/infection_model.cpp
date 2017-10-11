#include <Rcpp.h>
using namespace Rcpp;

#define MAX(a,b) ((a) < (b) ? (b) : (a)) // define MAX function for use later
  
// Do we want infection history to be a vector of infection times?
// We could have another vector that means we do not predict titres 
// when we are missing samples called:
// IntegerVector samples
//[[Rcpp::export]]
NumericVector infection_model_indiv(NumericVector theta, NumericVector infectionHistory,
				    double samplingTime, NumericVector strainIsolationTimes,
				    NumericVector antigenicMapLong,NumericVector antigenicMapShort){
  // Extract model parameters
  double mu = theta["mu"];
  double mu_short = theta["mu_short"];
  double tau = theta["tau"];
  double wane = theta["wane"];

  // We will need to loop over each strain that *could* be tested
  int n_strains = strainIsolationTimes.size();
  int n_samples = n_strains; // Number of time points sampled
  int max_infections = n_strains; // max number of infections is one for each strain
  //Rcpp::Rcout << "N strains indiv: " << n_strains << std::endl;
  double tmpTitre=0;

  NumericVector predictedTitre(n_samples);
  NumericVector cumInfectionHistory(max_infections);
  NumericVector maskedInfectionHistory(max_infections);
  NumericVector waning(n_samples);

  double circulation_time;

  // Set up cumulative infection history and waning vector
  /* Check if first isolation time is after the sampling time.
     if so, then we do not test against this strain */ 
  circulation_time = strainIsolationTimes[0];
  if(circulation_time > samplingTime) maskedInfectionHistory[0]=0;
  else maskedInfectionHistory[0] = infectionHistory[0];
  waning[0] = MAX(0, 1.0-wane*(samplingTime - circulation_time));
  cumInfectionHistory[0] = maskedInfectionHistory[0];

  /* For strains that circulated before the isolation time,
     add to cumulative infection history */
  //Rcpp::Rcout << "Sampling time: " << samplingTime << std::endl;
  for(int i=1; i < n_strains; ++i){
    /* At which time did this strain circulate?
       This might change to a range of times in the future */
    circulation_time = strainIsolationTimes[i];
    //Rcpp::Rcout << "Circulation time: " << circulation_time << std::endl;
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
    //Rcpp::Rcout << "Waning: " << waning[i] << std::endl;
  }
  // For each strain we are testing against, find predicted titre
  for(int k=0; k < n_samples; ++k){
    tmpTitre=0;
    /* Can add in the following code if we want to censor:
       if(samples[k] > 0){
       ...
       } else {
       predictedTitre[k] = NA_REAL;
       }
    */
    // Sum contributions from all infections
    for(int i=0; i < n_samples; ++i){
      ///////////////////////////////
      // THE ACTUAL MODEL
      tmpTitre += maskedInfectionHistory[i]* // Ignore infections that couldn't have happened
	MAX(0, 1.0 - tau*(cumInfectionHistory[i] - 1.0))* // Antigenic seniority
	(mu*antigenicMapLong[k*n_samples+i] +  // Long term boost
	 mu_short*antigenicMapShort[k*n_samples+i]* // Short term cross reactive boost
	 waning[i]); // Waning rate
      ////////////////////////////
    }
    predictedTitre[k] = tmpTitre;
  }
  return(predictedTitre);
}

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

//[[Rcpp::export]]
double individual_likelihood(NumericVector theta, NumericVector infectionHistory,
			     IntegerVector samples,
			     NumericVector samplingTimes, NumericVector strainIsolationTimes,
			     NumericVector antigenicMapLong, NumericVector antigenicMapShort, 
			     NumericVector titres){
  int numberStrains = strainIsolationTimes.size();
  int numberSamples = samplingTimes.size();
  NumericVector predictedTitres(numberStrains);

  double lnlike = 0;

  // These indices allow us to step through a vector as if it were a 
  // matrix (ie. numberStrains at a time)
  int startIndex = 0;
  int endIndex = numberStrains-1;
  
  for(int i=0; i < numberSamples; ++i){
    // Check if this individual actually has a sample for this time. If not, skip
    if(samples[i] > 0){
      predictedTitres = infection_model_indiv(theta, infectionHistory, samplingTimes[i],
					      strainIsolationTimes, antigenicMapLong, antigenicMapShort);
      lnlike += likelihood_titre(predictedTitres, titres[Range(startIndex,endIndex)],theta);
    }

    // Increase indices by max number of strains
    startIndex += numberStrains;
    endIndex += numberStrains;
  }
  return(lnlike);
}

//[[Rcpp::export]]
double group_likelihood(NumericVector theta, 
			NumericMatrix infectionHistories, IntegerVector samples,
			NumericVector samplingTimes, NumericVector strainIsolationTimes,
			NumericVector antigenicMapLong, NumericVector antigenicMapShort, 
			NumericVector titres){
  int n = infectionHistories.nrow();
  int n_strains = infectionHistories.ncol();
  int n_samples = samplingTimes.size();
  int indiv_length = n_strains*n_samples;

  IntegerVector tmpSamples(n_samples);
  
  // These indices allow us to step through the titre data vector
  // as if it were a matrix ie. number of rows for each individual
  // at a time
  int startIndex = 0;
  int endIndex = indiv_length-1;

  // There is a matrix that specifies if a sampling time is present
  // for each individual in each group. These indices allow us to step
  // through that matrix for each individual at a time
  int indivSampleIndexStart = 0;
  int indivSampleIndexEnd = n_samples-1;

  double lnlike = 0;

  for(int i=0; i < n; ++i){
    // Check if samples exist for this individual at each potential sampling time
    tmpSamples = samples[Range(indivSampleIndexStart, indivSampleIndexEnd)];
 
    lnlike += individual_likelihood(theta, infectionHistories(i,_), tmpSamples, samplingTimes, 
				    strainIsolationTimes, antigenicMapLong, antigenicMapShort, 
				    titres[Range(startIndex, endIndex)]);
    startIndex += indiv_length;
    endIndex += indiv_length;

    indivSampleIndexStart += n_samples;
    indivSampleIndexEnd += n_samples;
  }
  return(lnlike);  
}

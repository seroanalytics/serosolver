#include <Rcpp.h>
using namespace Rcpp;

#define MAX(a,b) ((a) < (b) ? (b) : (a)) // define MAX function for use later


//[[Rcpp::export]]
NumericVector subset_test(NumericVector x, IntegerVector y){
  LogicalVector indices = y > 0;
  return(x[indices]);
}


//[[Rcpp::export]]
NumericVector subset_test1(NumericVector x, LogicalVector y){
  return(x[y]);
}





//[[Rcpp::export]]
NumericVector infection_model_indiv_OLD(NumericVector theta, NumericVector infectionHistory,
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

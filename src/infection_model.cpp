#include <Rcpp.h>
#include <cmath>
#include "wane_function.h"
#include "boosting_functions.h"
using namespace Rcpp;

#define MAX(a,b) ((a) < (b) ? (b) : (a)) // define MAX function for use later


//' Marginal prior probability (p(Z)) of a particular infection history matrix
//' 
//' @param infHist IntegerMatrix, the infection history matrix
//' @param n_alive IntegerVector, vector giving the number of individuals alive in each time unit
//' @param alpha NumericVector, alpha parameter for beta distribution prior, one for each time unit
//' @param beta NumericVector, beta parameter for beta distribution prior, one for each time unit
//' @export
// [[Rcpp::export]]
double inf_mat_prior_cpp(const IntegerMatrix& infHist, const IntegerVector& n_alive, double alpha, double beta){
//double inf_mat_prior_cpp(const IntegerMatrix& infHist, const IntegerVector& n_alive, const NumericVector& alphas, const NumericVector& betas){
  // Prior on each year
  double m, n;
  double lik=0;
  double lbeta_const = R::lbeta(alpha, beta);
  for(int i = 0; i < n_alive.size(); ++i){
    m = sum(infHist(_,i));
    n = n_alive(i);
    lik += R::lbeta(m+alpha,n-m+beta)-lbeta_const;
  }
  return(lik);
}


//' Marginal prior probability (p(Z)) of a particular infection history matrix
//' 
//' @param infHist IntegerMatrix, the infection history matrix
//' @param n_alive IntegerVector, vector giving the number of individuals alive in each year
//' @param alpha double, alpha parameter for beta distribution prior
//' @param beta double, beta parameter for beta distribution prior
//' @export
// [[Rcpp::export]]
double inf_mat_prior_total_cpp(const IntegerMatrix& infHist, const int& n_alive, double alpha, double beta){
  double m, n;
  double lik=0;
  int n_infections = sum(infHist);
  lik = R::lbeta(n_infections + alpha, n_alive - n_infections + beta) - R::lbeta(alpha, beta);
  return(lik);
}


//' @export
// [[Rcpp::export(rng = false)]]
NumericVector likelihood_func_fast(NumericVector theta, NumericVector obs, NumericVector predicted_titres){
  int total_titres = predicted_titres.size();
  NumericVector ret(total_titres);
  const double sd = theta["error"];
  const double den = sd*M_SQRT2;
  const double max_titre = theta["MAX_TITRE"];
  const double log_const = log(0.5);

  for(int i = 0; i < total_titres; ++i){
    if(obs[i] <= max_titre && obs[i] > 0.0){
      ret[i] = log_const + log((erf((obs[i] + 1.0 - predicted_titres[i]) / den) -
				erf((obs[i]     - predicted_titres[i]) / den)));    
    } else if(obs[i] > max_titre) {
      ret[i] = log_const + log(1.0 + erf((max_titre - predicted_titres[i])/den));
    } else {
      ret[i] = log_const + log(1.0 + erf((1.0 - predicted_titres[i])/den));
    }
  }
  return(ret);
} 


double erf_native(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return sign*y;
}

//' @export
// [[Rcpp::export]]
NumericVector likelihood_func_fast_native(const NumericVector &theta, const NumericVector &obs, const NumericVector &predicted_titres){
  int total_titres = predicted_titres.size();
  NumericVector ret(total_titres);
  const double sd = theta["error"];
  const double den = sd*M_SQRT2;
  const double max_titre = theta["MAX_TITRE"];
  const double log_const = log(0.5);

  for(int i = 0; i < total_titres; ++i){
    if(obs[i] <= max_titre && obs[i] > 0.0){
      ret[i] = log_const + log((erf_native((obs[i] + 1.0 - predicted_titres[i]) / den) -
				erf_native((obs[i]     - predicted_titres[i]) / den)));    
    } else if(obs[i] > max_titre) {
      ret[i] = log_const + log(1.0 + erf_native((max_titre - predicted_titres[i])/den));
    } else {
      ret[i] = log_const + log(1.0 + erf_native((1.0 - predicted_titres[i])/den));
    }
  }
  return(ret);
} 





//' Calculate likelihood basic
//'
//' Calculates the likelihood of a given set of observed titres given predicted titres. Based on truncated discritised normal. DEPRECATED, as this is a very slow (but obvious) implementation
//' @param expected NumericVector, as returned by \code{\link{infection_model_indiv}}
//' @param data NumericVector, the vector of observed titres
//' @param theta NumericVector, the vector of named model parameters, requiring MAX_TITRE and error
//' @param titre_shifts NumericVector, OPTIONAL if using measurement bias, gives the shift to add to each expected titre
//' @return a single log likelihood
//' @export
//[[Rcpp::export]]
double likelihood_titre_basic(const NumericVector &expected, 
			const NumericVector &data, 
			const NumericVector &theta,
			const NumericVector &titre_shifts
			){

  NumericVector expected1 = expected;
  if(titre_shifts.size() == expected.size()){
    expected1 = expected1 + titre_shifts;
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
NumericVector sum_likelihoods(NumericVector liks, IntegerVector indices, int n_indivs){
  NumericVector results(n_indivs);
  int end = liks.size();
  for(int i = 0; i < end; ++i){
    results[indices[i]] += liks[i];
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

  
// Does some precomputation to speed up the model solving code
void setup_waning_and_masked_cumulative(const NumericVector &theta,
					const IntegerVector &infection_history,
					IntegerVector &cumu_infection_history, 
					IntegerVector &masked_infection_history,
					NumericVector &waning,
					NumericVector &seniority,
					const NumericVector &infection_times,
					const int &max_infections, 
					const double &sampling_time,
					const int &wane_type,
					const double DOB=0,
					const Nullable<List> &additional_arguments=R_NilValue){
  double circulation_time=infection_times[0];
  double time = sampling_time - circulation_time;
  double wane = theta["wane"];
  double tau = theta["tau"];
  double val;  
  
  /* =========================================
     ======== SET UP WANING RATE VECTOR =========
     ========================================= */ 
  // If not linear
  if(wane_type == 0){
    waning[0] = MAX(0, 1.0-wane*time);// Else if linear
  } else {
    // Calculate the waning at the time since infection
    val= wane_function(theta, time, wane);
    waning[0] = MAX(0, 1.0-val);
  }

  // Seniority
  seniority[0] = 1;

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
    if(wane_type == 0){
      waning[i] = MAX(0, 1.0-wane*time); // Else linear
    } else {
      val= wane_function(theta, time, wane);
      waning[i] = MAX(0, 1.0-val);
    }

    // Seniority
    seniority[i] = MAX(0, 1.0 - tau*(cumu_infection_history[i]-1));
  }  
}


//' Model function sample
//'
//' The main model solving function for a single individual for a single blood sample.
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
//' @param DOB double, the date of birth of this individual. Currently not used.
//' @param additional_arguments, Nullable<List> currently not used, but the idea is to use thsi object to pass more flexible additional arguments to the bottom of the call stack
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
  //IntegerVector cumInfectionHistory = cumsum(infectionHistory);

  // Trying just ignoring this, as subset in function above
  IntegerVector maskedInfectionHistory(max_infections); // To avoid solving the model for infections that didn't happen
  NumericVector waning(max_infections);
  NumericVector seniority(max_infections);
  
   // Set up cumulative infection history, masked infection history and waning vector
  setup_waning_and_masked_cumulative(theta,
				     infectionHistory,
				     cumInfectionHistory, 
				     maskedInfectionHistory, 
				     waning,	
				     seniority,			     
				     infectionTimes, 
				     max_infections, 
				     samplingTime,
				     waneType, 
				     DOB,
				     additional_arguments);

  add_multiple_infections_boost(predictedTitre, monitored_titres,
				theta, 
				infectionTimes, cumInfectionHistory, 
				maskedInfectionHistory, 
				infectionMapIndices, measurementMapIndices, 
				antigenicMapLong, antigenicMapShort,
				waning, 
				seniority,
				numberStrains, 
				n_samples,
				max_infections,
				titre_dependent_boosting, 
				DOB,
				additional_arguments
				);
  return(predictedTitre);
}

//' Model function individual
//'
//' The main model solving function for a single individual for a vector of sampling times
//' @param theta NumericVector, the named vector of model parameters
//' @param infectionHistory IntegerVector, the vector of 1s and 0s showing presence/absence of infection for each possible time. 
//' @param circulationTimes NumericVector, the actual times of circulation that the infection history vector corresponds to
//' @param circualtionaMapIndices IntegerVector, which entry in the melted antigenic map that these infection times correspond to
//' @param samplingTime NumericVector, the times that each blood sample was taken
//' @param dataIndices IntegerVector, the indices in the overall titre data vector (of observations) that each sample corresponds to
//' @param measuredMapIndices IntegerVector, the indices of all measured strains in the melted antigenic map
//' @param antigenicMapLong NumericVector, the collapsed cross reactivity map for long term boosting, after multiplying by sigma1
//' @param antigenicMapShort NumericVector, the collapsed cross reactivity map for short term boosting, after multiplying by sigma2
//' @param numberStrains int, the maximum number of infections that an individual could experience
//' @param DOB double, the date of birth of this individual. Currently not used.
//' @param additional_arguments, Nullable<List> currently not used, but the idea is to use thsi object to pass more flexible additional arguments to the bottom of the call stack
//' @return NumericVector of predicted titres for each entry in measuredMapIndices
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

  IntegerVector tmp_range;

  for(int i = 0; i < numberSamples; ++i){
    endIndex = startIndex + dataIndices[i] - 1;
    // Range index twice
    tmp_range = Range(startIndex, endIndex);
    titres[tmp_range] = infection_model_indiv(theta,conciseInfHist,infectionTimes,infMapIndices,
					       samplingTimes[i],measuredMapIndices[tmp_range],
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
  IntegerVector tmp_range;
  IntegerVector tmp_range_samples;
  int DOB;

  for(int i=1; i <= n; ++i){
    startIndexSamples = indicesSamples[i-1];
    endIndexSamples = indicesSamples[i] - 1;

    startIndexData = indicesTitreDataOverall[i-1];
    endIndexData = indicesTitreDataOverall[i] - 1;

    DOB = DOBs[i-1];

    tmp_range = Range(startIndexData, endIndexData);
    tmp_range_samples = Range(startIndexSamples, endIndexSamples);

    titres[tmp_range] = titre_data_individual(theta,   // Vector of named model parameters
					      infectionHistories(i-1,_), // Vector of infection history for individual i
					      circulationTimes, // Vector of all virus circulation times, same length and ncol infectionHistories
					      circulationMapIndices, // Gives the corresponding index in the antigenic map vector
					      samplingTimes[tmp_range_samples],  // Get sampling times for this individual
					      indicesTitreDataSample[tmp_range_samples], // The 
					      measuredMapIndices[tmp_range],  // For the indices in the antigenic map to which each titre corresponds
					      antigenicMapLong, 
					      antigenicMapShort,  
					      n_strains,
					      DOB, 
					      additional_arguments); // The total number of strains that circulated
    
  }
  return(titres);
}


void titre_data_fast_individual_base(NumericVector &predicted_titres,
				     const double &mu, const double &mu_short, 
				     const double &wane, const double &tau,
				     const NumericVector &infection_times,
				     const IntegerVector &infection_strain_indices_tmp,
				     const IntegerVector &measurement_strain_indices,
				     const NumericVector &sample_times,
				     const int &index_in_samples,
				     const int &end_index_in_samples,
				     const int &start_index_in_data,
				     const IntegerVector &nrows_per_blood_sample,
				     const int &number_strains,
				     const NumericVector &antigenic_map_short,
				     const NumericVector &antigenic_map_long
				     ){
  double sampling_time;
  double time;
  double n_inf;
  double wane_amount;
  double seniority;
  
  int n_titres;
  int max_infections = infection_times.size();
  int end_index_in_data;
  int tmp_titre_index = start_index_in_data;
  int inf_map_index;
  int index;    

  for(int j = index_in_samples; j <= end_index_in_samples; ++j){
    sampling_time = sample_times[j];
    n_inf = 1.0;	
    n_titres = nrows_per_blood_sample[j];
	
    end_index_in_data = start_index_in_data + n_titres;
    tmp_titre_index = start_index_in_data;

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
    tmp_titre_index = end_index_in_data;
  }

}
				     


//' @export
// [[Rcpp::export(rng = false)]]
NumericVector titre_data_fast(const NumericVector &theta, 
			      const IntegerMatrix &infection_history_mat, 
			      const NumericVector &circulation_times,
			      const IntegerVector &circulation_times_indices,
			      const NumericVector &sample_times,
			      const IntegerVector &rows_per_indiv_in_samples, // How many rows in titre data correspond to each individual, sample and repeat?
			      const IntegerVector &cum_nrows_per_individual_in_data, // How many rows in the titre data correspond to each individual?
			      const IntegerVector &nrows_per_blood_sample, // Split the sample times and runs for each individual
			      const IntegerVector &measurement_strain_indices, // For each titre measurement, corresponding entry in antigenic map
			      const NumericVector &antigenic_map_long,
			      const NumericVector &antigenic_map_short
			      ){
  int n = infection_history_mat.nrow();
  int number_strains = infection_history_mat.ncol();
  int total_titres = measurement_strain_indices.size();
  int max_infections;
  int n_titres;

  int index_in_samples;
  int end_index_in_samples;
  int number_samples;
  int start_index_in_data;
  int end_index_in_data;
  int tmp_titre_index;
  int inf_map_index;
  int index;
  double sampling_time;
  double time;
  double n_inf;

  IntegerVector infection_history(number_strains);
  LogicalVector indices;

  NumericVector infection_times;
  IntegerVector infection_strain_indices_tmp;

  double mu = theta["mu"];
  double mu_short = theta["mu_short"];
  double wane = theta["wane"];
  double tau = theta["tau"];
  double wane_amount;
  double seniority;

  NumericVector predicted_titres(total_titres);

  for(int i = 1; i <= n; ++i){
    infection_history = infection_history_mat(i-1,_);
    indices = infection_history > 0;
    infection_times = circulation_times[indices];

    if(infection_times.size() > 0){
      infection_strain_indices_tmp = circulation_times_indices[indices];
      max_infections = infection_times.size();

      index_in_samples = rows_per_indiv_in_samples[i-1];
      end_index_in_samples = rows_per_indiv_in_samples[i] - 1;
      number_samples = end_index_in_samples - index_in_samples;      
      start_index_in_data = cum_nrows_per_individual_in_data[i-1];

      for(int j = index_in_samples; j <= end_index_in_samples; ++j){
	sampling_time = sample_times[j];
	n_inf = 1.0;	
	n_titres = nrows_per_blood_sample[j];
	
	end_index_in_data = start_index_in_data + n_titres;
	tmp_titre_index = start_index_in_data;

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
  }  
  return(predicted_titres);
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
				  const NumericVector &titre_shifts,
				  const double &DOB,
				  const Nullable<List> &additional_arguments
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
  lnlike = likelihood_titre_basic(titres, data, theta, titre_shifts);
  return(lnlike);
}


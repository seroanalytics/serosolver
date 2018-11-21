#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
#define MAX(a,b) ((a) < (b) ? (b) : (a)) // define MAX function for use later

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

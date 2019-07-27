#include "helpers.h"

//' Takes a subset of a Nullable NumericVector, but only if it isn't NULL
// [[Rcpp::export]]
NumericVector subset_nullable_vector(const Nullable<NumericVector> &x, int index1, int index2) {
  if(x.isNotNull()){
    NumericVector y = as<NumericVector>(x)[Range(index1, index2)];
    return y;
  } else {
    NumericVector y(0);
    return y;
  }
}

//' Sum likelihoods into buckets
//' 
//' Given a large vector of likelihood values and a vector, indices, of length n_indivs, sums the likelihoods to give one value per individual, where indices indicates which individual each index of liks corresponds to.
//' @param liks NumericVector of likelihoods
//' @param indices IntegerVector of indices of same length as liks, where the max value of this should be the same as n_indivs - 1
//' @param n_indivs int, number of individuals to generate bucketed likelihoods for
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

//' Count infections by group and time
//'
//' @export
//[[Rcpp::export]]
IntegerMatrix sum_infections_by_group(IntegerMatrix inf_hist, IntegerVector group_ids_vec, int n_groups){
  int n_times = inf_hist.ncol();
  IntegerMatrix n_infections(n_groups, n_times);

  for(int i = 0; i < n_times; ++i){
    for(int j = 0; j < group_ids_vec.size(); ++j){
      n_infections(group_ids_vec[j], i) += inf_hist(j, i);
    }   
  }
  return(n_infections);
}

//' Add measurement shifts to predictions
//'
//' Adds observation error shifts to predicted antibody titres.
//' @param predicted_titre NumericVector, the predicted titres. Note that this vector will be changed!
//' @param to_add NumericVector the vector of all measurement shifts to apply
//' @param start_index_in_data int the first index of to_add and predicted_titres to combine
//' @param end_index_in_data int the end index of to_add and predicted_titres to combine
//' @return nothing
//' @export
//[[Rcpp::export]]
void add_measurement_shifts(NumericVector &predicted_titres, 
			    const NumericVector &to_add,
			    const int &start_index_in_data,
			    const int &end_index_in_data
			    ){
  for(int j = start_index_in_data; j <= end_index_in_data; ++j){
    predicted_titres[j] += to_add[j];
  }
}
			    



//' Titre protection
//'
//' @export
//[[Rcpp::export]]
double titre_protection_cpp(double titre, double alpha1, double beta1){
  double risk = 1.0 - 1.0/(1.0 + exp(beta1*(titre - alpha1)));
  return risk;
}

//' Titre protection
//'
//' @export
//[[Rcpp::export]]
double p_infection_cpp(double phi, double titre, double alpha1, double beta1){
  double p = phi*(1.0 - titre_protection_cpp(titre, alpha1, beta1));
  return p;
}



//' Calc titre probs
//'
//' @export
//[[Rcpp::export]]
NumericVector calc_phi_probs_indiv_titre_cpp(NumericVector phis, NumericVector titres,
					     IntegerMatrix infection_history,
					     IntegerVector age_mask,
					     IntegerVector strain_mask,
					     double alpha1, double beta1){
  int n_indivs = age_mask.size();
  int n_strains = infection_history.ncol();
  int index = 0;
  int x;
  double titre_p;
  double prob;
  double first;
  double second;
  NumericVector liks(n_indivs);
  for(int i = 0; i < n_indivs; ++i){
    for(int j = 0; j < n_strains; ++j){
      x = infection_history(i,j);
      index = i*n_strains + j;
      prob = 0;
      if(age_mask[i] <= (j+1) && strain_mask[i] >= (j+1)){
	titre_p = titre_protection_cpp(titres[index], alpha1, beta1);
	prob = log(x*(1-titre_p)*phis[j] + (1-x)*(1 + phis[j]*(titre_p*phis[j] - 1)));
      }
      liks[i] += prob;
    }
  }
  return liks;}


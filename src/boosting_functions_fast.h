#include <Rcpp.h>
using namespace Rcpp;

#ifndef TITRE_DATA_FAST_INDIVIDUAL_BASE_H
#define TITRE_DATA_FAST_INDIVIDUAL_BASE_H
void titre_data_fast_individual_base(NumericVector &predicted_titres,
				     const double &mu, const double &mu_short, 
				     const double &wane, const double &tau,
				     const NumericVector &wane_amounts,
				     const NumericVector &seniority_amounts,
				     const NumericVector &infection_times,
				     const IntegerVector &infection_strain_indices_tmp,
				     const IntegerVector &measurement_strain_indices,
				     const NumericVector &sample_times,
				     const int &index_in_samples,
				     const int &end_index_in_samples,
				     const int &start_index_in_data1,
				     const IntegerVector &nrows_per_blood_sample,
				     const int &number_strains,
				     const NumericVector &antigenic_map_short,
				     const NumericVector &antigenic_map_long,
				     bool boost_before_infection
				     );
#endif


#ifndef TITRE_DATA_FAST_INDIVIDUAL_TITREDEP_H
#define TITRE_DATA_FAST_INDIVIDUAL_TITREDEP_H
void titre_data_fast_individual_titredep(NumericVector &predicted_titres,
					   const double &mu,
					   const double &mu_short,
					   const double &wane,
					   const double &tau,
					   const double &gradient,
					   const double &boost_limit,
					   const NumericVector &infection_times,
					   const IntegerVector &infection_strain_indices_tmp,
					   const IntegerVector &measurement_strain_indices,
					   const NumericVector &sample_times,
					   const int &index_in_samples,
					   const int &end_index_in_samples,
					   const int &start_index_in_data1,
					   const IntegerVector &nrows_per_blood_sample,
					   const int &number_strains,
					   const NumericVector &antigenic_map_short,
					 const NumericVector &antigenic_map_long,
					 bool boost_before_infection
					 );
#endif


#ifndef TITRE_DATA_FAST_INDIVIDUAL_COMPLEX_CR_H
#define TITRE_DATA_FAST_INDIVIDUAL_COMPLEX_CR_H
void titre_data_fast_individual_complex_cr(NumericVector &predicted_titres,
                                           const double &mu,
                                           const double &mu_short,
                                           const double &wane,
                                           const double &tau,
                                           
                                           // Cross reactivity parameters
                                           const double &sigma_s, // short term cross reactivity
                                           const double &sigma_l, // long term cross reactivity
                                           const double &sigma_birth_mod_s, // short term cross-reactivity modifier to pre-birth
                                           const double &sigma_birth_mod_l, // long term cross-reactivity modifier to pre-birth
                                           const double &sigma_future_mod_s, // short term cross-reactivity modifier to future
                                           const double &sigma_future_mod_l, // long term cross-reactivity modifier to future
                                           
                                           const NumericVector &wane_amounts,
                                           const NumericVector &seniority_amounts,
                                           const NumericVector &infection_times,
                                           const IntegerVector &infection_strain_indices_tmp,
                                           const IntegerVector &measurement_strain_indices,
                                           const NumericVector &sample_times,
                                           const int &index_in_samples,
                                           const int &end_index_in_samples,
                                           const int &start_index_in_data1,
                                           const IntegerVector &nrows_per_blood_sample,
                                           const int &number_strains,
                                           const NumericVector &antigenic_distances,
                                           bool boost_before_infection
);
#endif
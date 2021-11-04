#include <Rcpp.h>
using namespace Rcpp;

#ifndef TITRE_DATA_FAST_INDIVIDUAL_BASE_H
#define TITRE_DATA_FAST_INDIVIDUAL_BASE_H
void titre_data_fast_individual_base(NumericVector &predicted_titres,
				     const double &mu, const double &mu_short, 
				     const double &wane, const double &tau,
                                     const bool &vac_flag,
                                    // const IntegerVector &vac_flag_ind,
                                     const double &mu_vac,
                                     const double &mu_short_vac,
                                     const double &wane_vac,
									 const double &tau_prev_vac,
                     NumericVector infection_times,
                     IntegerVector infection_strain_indices_tmp,
                     const NumericVector &vaccination_times,
                     const IntegerVector &vaccination_strain_indices_tmp,
					 const NumericVector &vaccination_previous,
				     const IntegerVector &measurement_strain_indices,
				     const NumericVector &sample_times,
				     const int &index_in_samples,
				     const int &end_index_in_samples,
				     const int &start_index_in_data1,
				     const IntegerVector &nrows_per_blood_sample,
				     const int &number_strains,
				     const NumericVector &antigenic_map_short,
				     const NumericVector &antigenic_map_long,
					 const NumericVector &antigenic_map_short_vac,
				     const NumericVector &antigenic_map_long_vac,
				     bool boost_before_infection
				     );
#endif


#ifndef TITRE_DATA_FAST_INDIVIDUAL_WANE2_H
#define TITRE_DATA_FAST_INDIVIDUAL_WANE2_H
void titre_data_fast_individual_wane2(NumericVector &predicted_titres,
				      const double &mu,
				      const double &mu_short,
				      const double &wane,
				      const double &tau,
				      const double &kappa,
				      const double &t_change,
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

#ifndef TITRE_DATA_FAST_INDIVIDUAL_STRAIN_DEPENDENT_H
#define TITRE_DATA_FAST_INDIVIDUAL_STRAIN_DEPENDENT_H
void titre_data_fast_individual_strain_dependent(NumericVector &predicted_titres,
						 const NumericVector &mus,
						 const IntegerVector &boosting_vec_indices,
						 const double &mu_short,
						 const double &wane,
						 const double &tau,
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

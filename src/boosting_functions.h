#include <Rcpp.h>
using namespace Rcpp;


#ifndef MULTIPLE_INFECTION_BASE_BOOSTING_H
#define MULTIPLE_INFECTION_BASE_BOOSTING_H
void multiple_infection_base_boosting(NumericVector &predicted_titres,
				      const NumericVector &theta,
				      const IntegerVector &cumu_infection_history,
				      const IntegerVector &masked_infection_history,
				      const IntegerVector &infection_map_indices, 
				      const IntegerVector &measurement_map_indices,
				      const NumericVector &antigenic_map_long, 
				      const NumericVector &antigenic_map_short, 
				      const NumericVector &waning,
				      const NumericVector &seniority,
				      const int &number_strains,
				      const int &n_samples,
				      const int &max_infections
				      );

#endif

#ifndef MULTIPLE_INFECTION_TITRE_DEPENDENT_BOOST_H
#define MULTIPLE_INFECTION_TITRE_DEPENDENT_BOOST_H
void multiple_infection_titre_dependent_boost(NumericVector &predicted_titres, 
					      NumericVector &monitored_titres,
					      const NumericVector &theta,
					      const NumericVector &infection_times,
					      const IntegerVector &cumu_infection_history,
					      const IntegerVector &masked_infection_history,
					      const IntegerVector &infection_map_indices,
					      const IntegerVector &measurement_map_indices,
					      const NumericVector &antigenic_map_long, 
					      const NumericVector &antigenic_map_short, 
					      const NumericVector &waning,
					      const int &number_strains
					      );
#endif

#ifndef MULTIPLE_INFECTION_STRAIN_DEPENDENT_H
#define MULTIPLE_INFECTION_STRAIN_DEPENDENT_H
void multiple_infection_strain_dependent(NumericVector &predicted_titres,
					 const NumericVector &theta,
					 const IntegerVector &cumu_infection_history,
					 const IntegerVector &masked_infection_history,
					 const IntegerVector &infection_map_indices, 
					 const IntegerVector &measurement_map_indices,
					 const NumericVector &antigenic_map_long, 
					 const NumericVector &antigenic_map_short, 
					 const NumericVector &waning,
					 const int &number_strains,
					 List additional_arguments);
#endif


#ifndef ADD_MULTIPLE_INFECTIONS_BOOST_H
#define ADD_MULTIPLE_INFECTIONS_BOOST_H
void add_multiple_infections_boost(NumericVector &predicted_titres, 
				   NumericVector &monitored_titres,
				   const NumericVector &theta,
				   const NumericVector &infection_times,
				   const IntegerVector &cumu_infection_history,
				   const IntegerVector &masked_infection_history,
				   const IntegerVector &infection_map_indices,
				   const IntegerVector &measurement_map_indices,
				   const NumericVector &antigenic_map_long, 
				   const NumericVector &antigenic_map_short, 
				   const NumericVector &waning,
				   const NumericVector &seniority,
				   const int &number_strains,
				   const int &n_samples,
				   const int &max_infections,
				   const bool &titre_dependent_boosting,
				   const int &DOB,
				   const Nullable<List> &additional_arguments
				   );
#endif

#ifndef TITRE_DATA_FAST_INDIVIDUAL_BASE_H
#define TITRE_DATA_FAST_INDIVIDUAL_BASE_H
void titre_data_fast_individual_base(NumericVector &predicted_titres,
				     const double &mu, const double &mu_short, 
				     const double &wane, const double &tau,
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
				     const NumericVector &antigenic_map_long
				     );
#endif

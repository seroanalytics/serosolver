#include <Rcpp.h>
using namespace Rcpp;


#ifndef MULTIPLE_INFECTION_BASE_BOOSTING_H
#define MULTIPLE_INFECTION_BASE_BOOSTING_H
void multiple_infection_base_boosting(NumericVector &predicted_titres,
				      const NumericVector &theta,
				      const IntegerVector &cumu_infection_history,
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

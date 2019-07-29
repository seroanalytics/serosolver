#include <Rcpp.h>
using namespace Rcpp;

#ifndef LIKELIHOOD_TITRE_BASIC_H
#define LIKELIHOOD_TITRE_BASIC_H
double likelihood_titre_basic(const NumericVector &expected, 
			      const NumericVector &data, 
			      const NumericVector &theta,
			      const NumericVector &titre_shifts
			      );
#endif

#ifndef PROPOSAL_LIKELIHOOD_FUNC_H
#define PROPOSAL_LIKELIHOOD_FUNC_H
void proposal_likelihood_func(double &new_prob,
			      NumericVector &predicted_titres,
			      const int &indiv,
			      const NumericVector &data,
			      const NumericVector &repeat_data,
			      const IntegerVector &repeat_indices,
			      const IntegerVector &cum_nrows_per_individual_in_data,
			      const IntegerVector &cum_nrows_per_individual_in_repeat_data,
			      const double &log_const,
			      const double &den,
			      const double &max_titre,
			      const bool &repeat_data_exist);
#endif

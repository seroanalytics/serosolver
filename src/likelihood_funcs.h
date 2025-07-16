#include <Rcpp.h>
using namespace Rcpp;

#ifndef PROPOSAL_LIKELIHOOD_FUNC_H
#define PROPOSAL_LIKELIHOOD_FUNC_H
void proposal_likelihood_func(double &new_prob,
                              NumericVector &predicted_antibody_levels,
                              const int &indiv,
                              const NumericVector &data,
                              const NumericVector &repeat_data,
                              const IntegerVector &repeat_indices,
                              const IntegerVector &cum_nrows_per_individual_in_data,
                              const IntegerVector &cum_nrows_per_individual_in_repeat_data,
                              const double &log_const,
                              const double &den,
                              const double &max_measurement,
                              const double &min_measurement,
                              const bool &repeat_data_exist,
                              const double &obs_weight);
#endif

#ifndef PROPOSAL_LIKELIHOOD_FUNC_CONTINUOUS_H
#define PROPOSAL_LIKELIHOOD_FUNC_CONTINUOUS_H
void proposal_likelihood_func_continuous(double &new_prob,
                                         NumericVector &predicted_antibody_levels,
                                         const int &indiv,
                                         const NumericVector &data,
                                         const NumericVector &repeat_data,
                                         const IntegerVector &repeat_indices,
                                         const IntegerVector &cum_nrows_per_individual_in_data,
                                         const IntegerVector &cum_nrows_per_individual_in_repeat_data,
                                         const double &log_const,
                                         const double &sd,
                                         const double &den,
                                         const double &den2,
                                         const double &max_measurement,
                                         const double &min_measurement,
                                         const bool &repeat_data_exist,
                                         const double &obs_weight);
#endif

#ifndef PROPOSAL_LIKELIHOOD_FUNC_CONTINUOUS_FP_H
#define PROPOSAL_LIKELIHOOD_FUNC_CONTINUOUS_FP_H
void proposal_likelihood_func_continuous_fp(double &new_prob,
                                         NumericVector &predicted_antibody_levels,
                                         const int &indiv,
                                         const int &prev_infected_flag,
                                         const NumericVector &data,
                                         const NumericVector &repeat_data,
                                         const IntegerVector &repeat_indices,
                                         const IntegerVector &cum_nrows_per_individual_in_data,
                                         const IntegerVector &cum_nrows_per_individual_in_repeat_data,
                                         const double &log_const,
                                         const double &sd,
                                         const double &den,
                                         const double &den2,
                                         const double &log_fp_rate,
                                         const double &log_fp_rate_compliment,
                                         const double &log_unif_pdf,
                                         const double &max_measurement,
                                         const double &min_measurement,
                                         const bool &repeat_data_exist,
                                         const double &obs_weight);
#endif


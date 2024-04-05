#include <RcppArmadillo.h>
using namespace Rcpp;

#ifndef ANTIBODY_DEPENDENT_BOOSTING_MODEL_INDIVIDUAL_H
#define ANTIBODY_DEPENDENT_BOOSTING_MODEL_INDIVIDUAL_H
void antibody_dependent_boosting_model_individual(NumericVector &predicted_antibody_levels,
                                                  const double &boost_long,
                                                  const double &boost_short,
                                                  const double &wane_short,
                                                  const double &antigenic_seniority,
                                                  const double &gradient,
                                                  const double &boost_limit,
                                                  const NumericVector &possible_exposure_times,
                                                  const IntegerVector &exposure_indices,
                                                  const IntegerVector &biomarker_id_indices,
                                                  const NumericVector &sample_times,
                                                  const int &index_in_samples,
                                                  const int &end_index_in_samples,
                                                  const int &start_index_in_data1,
                                                  const IntegerVector &nrows_per_blood_sample,
                                                  const int &number_possible_exposures,
                                                  const double *antigenic_map_short,
                                                  const double *antigenic_map_long,
                                                  bool boost_before_infection
);
#endif


#ifndef ANTIBODY_DATA_MODEL_INDIVIDUAL_NEW_H
#define ANTIBODY_DATA_MODEL_INDIVIDUAL_NEW_H
void antibody_data_model_individual_new(NumericVector &predicted_antibody_levels,
                                        const NumericVector &start_antibody_levels,
                                        const NumericVector &births,
                                        const double &boost_long,
                                        const double &boost_short,
                                        const double &boost_delay,
                                        const double &wane_short,
                                        const double &wane_long,
                                        const double &antigenic_seniority,
                                        const NumericVector &possible_exposure_times,
                                        const IntegerVector &exposure_indices,
                                        const IntegerVector &biomarker_id_indices,
                                        const IntegerVector &start_level_indices,
                                        const NumericVector &sample_times,
                                        const int &index_in_samples,
                                        const int &end_index_in_samples,
                                        const int &start_index_in_data1,
                                        const IntegerVector &nrows_per_blood_sample,
                                        const int &number_possible_exposures,
                                        const double *antigenic_map_short,
                                        const double *antigenic_map_long,
                                        bool boost_before_infection,
                                        const double min_level
);
#endif

#ifndef ANTIBODY_DATA_MODEL_INDIVIDUAL_H
#define ANTIBODY_DATA_MODEL_INDIVIDUAL_H
void antibody_data_model_individual(NumericVector &predicted_antibody_levels,
                                    const double &boost_long,
                                    const double &boost_short,
                                    const double &boost_delay,
                                    const double &wane_short,
                                    const double &wane_long,
                                    const double &antigenic_seniority,
                                    const NumericVector &possible_exposure_times,
                                    const IntegerVector &exposure_indices,
                                    const IntegerVector &biomarker_id_indices,
                                    const NumericVector &sample_times,
                                    const int &index_in_samples,
                                    const int &end_index_in_samples,
                                    const int &start_index_in_data1,
                                    const IntegerVector &nrows_per_blood_sample,
                                    const int &number_possible_exposures,
                                    const double *antigenic_map_short,
                                    const double *antigenic_map_long,
                                    bool boost_before_infection
);
#endif



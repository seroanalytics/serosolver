#include <RcppArmadillo.h>
using namespace Rcpp;

#ifndef MAX
#define MAX(a,b) ((a) < (b) ? (b) : (a)) // define MAX function for use later
#endif

#ifndef SUBSET_NULLABLE_VECTOR_H
#define SUBSET_NULLABLE_VECTOR_H
NumericVector subset_nullable_vector(const Nullable<NumericVector> &x, int index1, int index2);
#endif

#ifndef GET_STARTING_ANTIBODY_LEVELS_H
#define GET_STARTING_ANTIBODY_LEVELS_H
NumericVector get_starting_antibody_levels(const int n_measurements, 
                                           const double min_measurement, 
                                           const Nullable<NumericVector> &starting_antibody_levels);
#endif
  
#ifndef CREATE_CROSS_REACTIVITY_VECTOR_H
#define SUBSET_NULLABLE_VECTOR_H
NumericVector create_cross_reactivity_vector(NumericVector x, double cr_gradient);
#endif

#ifndef SUM_LIKELIHOODS_H
#define SUM_LIKELIHOODS_VECTOR_H
NumericVector sum_likelihoods(NumericVector liks, IntegerVector indices, int n_indivs);
#endif

#ifndef SUM_BUCKETS_H
#define SUM_BUCKETS_H
NumericVector sum_buckets(NumericVector a, NumericVector buckets);
#endif

#ifndef ADD_MEASUREMENT_SHIFTS
#define ADD_MEASUREMENT_SHIFTS
void add_measurement_shifts(NumericVector &predicted_titres, 
			    const NumericVector &to_add,
			    const int &start_index_in_data,
			    const int &end_index_in_data
			    );
#endif

#include <Rcpp.h>
using namespace Rcpp;


#ifndef MAX
#define MAX(a,b) ((a) < (b) ? (b) : (a)) // define MAX function for use later
#endif

#ifndef SUBSET_NULLABLE_VECTOR_H
#define SUBSET_NULLABLE_VECTOR_H
NumericVector subset_nullable_vector(const Nullable<NumericVector> &x, int index1, int index2);
#endif

#ifndef CREATE_CROSS_REACTIVITY_VECTOR_H
#define SUBSET_NULLABLE_VECTOR_H
NumericVector create_cross_reactivity_vector(NumericVector x, double sigma);
#endif

#ifndef SUM_LIKELIHOODS_H
#define SUM_LIKELIHOODS_VECTOR_H
NumericVector sum_likelihoods(NumericVector liks, IntegerVector indices, int n_indivs);
#endif

#ifndef SUM_BUCKETS_H
#define SUM_BUCKETS_H
NumericVector sum_buckets(NumericVector a, NumericVector buckets);
#endif

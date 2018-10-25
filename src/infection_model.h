#include <Rcpp.h>
using namespace Rcpp;


#ifndef LIKELIHOOD_DATA_INDIVIDUAL_H
#define LIKELIHOOD_DATA_INDIVIDUAL_H
double likelihood_data_individual(const NumericVector &theta, 
				  const IntegerVector &infectionHistory, 
				  const NumericVector &circulationTimes, 
				  const IntegerVector &circulationMapIndices,
				  const NumericVector &samplingTimes,
				  const IntegerVector &dataIndices,
				  const IntegerVector &measuredMapIndices, 
				  const NumericVector &antigenicMapLong, 
				  const NumericVector &antigenicMapShort,
				  const int &numberStrains,
				  const NumericVector &data,
				  const NumericVector &to_add,
				  const double &age,
				  const Nullable<List> &additional_arguments
				  ) ;
#endif

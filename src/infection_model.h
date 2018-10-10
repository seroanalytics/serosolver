#ifndef LIKELIHOOD_DATA_INDIVIDUAL_H
#define LIKELIHOOD_DATA_INDIVIDUAL_H


#include <Rcpp.h>
using namespace Rcpp;

double likelihood_data_individual(NumericVector theta, 
				  IntegerVector infectionHistory, 
				  NumericVector circulationTimes, 
				  IntegerVector circulationMapIndices,
				  NumericVector samplingTimes,
				  IntegerVector dataIndices,
				  IntegerVector measuredMapIndices, 
				  NumericVector antigenicMapLong, 
				  NumericVector antigenicMapShort,
				  int numberStrains,
				  NumericVector data,
				  Nullable<NumericVector> to_add
				  ) ;

#endif

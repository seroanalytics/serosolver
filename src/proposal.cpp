#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

#define MAX(a,b) ((a) < (b) ? (b) : (a)) // define MAX function for use later

//' Fast infection history proposal function
//' 
//' Proposes a new matrix of infection histories using a beta binomial proposal distribution. This particular implementation allows for nInfs epoch times to be changed with each function call. Furthermore, the size of the swap step is specified for each individual by moveSizes.
//' @param infHist and RcppArmadillo matrix of infection histories, where rows represent individuals and columns represent potential infection times. The contents should be a set of 1s (presence of infection) and 0s (absence of infection)
//' @param sampledIndivs IntegerVector, indices of which individuals to resample. Note that this is indexed from 1 (ie. as if passing straight from R)
//' @param ageMask IntegerVector, for each individual gives the first column in the infection history matrix that an individual could have been exposed to indexed from 1. ie. if alive for the whole period, entry would be 1. If alive for the 11th epoch, entry would be 11.
//' @param moveSizes IntegerVector, how far can a swap step sample from specified for each individual
//' @param nInfs IntegetVector, how many infections to add/remove/swap with each proposal step for each individual
//' @param alpha double, alpha parameter of the beta binomial
//' @param beta double, beta parameter of the beta binomial
//' @param randNs NumericVector, a vector of random numbers for each sampled individual. The idea is to pre-specify whether an individual experiences an add/remove step or a swap step to avoid random number sampling in C++
//' @return a matrix of 1s and 0s corresponding to the infection histories for all individuals
// [[Rcpp::export]]
arma::mat inf_hist_prop_cpp(arma::mat infHist, 
			    IntegerVector sampledIndivs, 
			    IntegerVector ageMask,
			    IntegerVector moveSizes, 
			    IntegerVector nInfs,
			    double alpha, 
			    double beta, 
			    NumericVector randNs) {
  // Copy input matrix
  arma::mat newInfHist = infHist;
  IntegerVector locs; // Locations to be updated
  arma::uvec locs1;
  arma::mat x;
  arma::mat y;
  IntegerVector samps;
  
  int index = 0;
  
  int maxI_indiv;
  int maxI = newInfHist.n_cols;
  int indiv;
  int k;
  int nInf;
  int n;
  int moveMax;
  int move;
  int id1;
  int id2;
  int tmp;
  
  double rand1;
  double ratio;
  
  // For each sampled individual
  for(int i = 0; i < sampledIndivs.size(); ++i){

    // Isolate that individual's infection histories
    indiv = sampledIndivs[i]-1;
    //Rcpp::Rcout << "individual: " << indiv << std::endl;
    nInf = nInfs[indiv];
    x = newInfHist.submat(indiv, ageMask[indiv]-1, indiv, maxI-1);
    samps = seq_len(x.n_cols);
    // With 50% probability, add/remove infections or swap infections
    if(randNs[i] < 1.0/2.0){
      /*
      Rcpp::Rcout << "Indiv: " << indiv << std::endl;
      Rcpp::Rcout << "ageMask: " << ageMask[indiv]-1 << std::endl;
      Rcpp::Rcout << "nInf: " << nInf << std::endl;
      Rcpp::Rcout << "Inf hist: " << x << std::endl;
      Rcpp::Rcout << "Samps: " << samps << std::endl;
      */
      // Sample N random locations
      locs = RcppArmadillo::sample(samps, nInf, FALSE, NumericVector::create());
      //Rcpp::Rcout << "Available locs: " << locs << std::endl;
      locs1 = as<arma::uvec>(locs)-1;
      //Rcpp::Rcout << "Locs chosen: " << locs1 << std::endl;
      //y = newInfHist.row(indiv);
      //Rcpp::Rcout << "y: " << y << std::endl;
      y = x.elem(locs1);
      //Rcpp::Rcout << "locations chosen contents: " << y << std::endl;
      // Count the number of 1s and 0s
      k = accu(x) - accu(y);
      //Rcpp::Rcout << "Infections less sampled: " << k << std::endl;
      n = x.size() - nInf;
      //Rcpp::Rcout << "N less sampled: " << n << std::endl;

      
      // For each sampled location, choose to turn into a 1 or 0 depending
      // on the beta binomial distribution.
      for(int j = 0; j < nInf; ++j){
	//Rcpp::Rcout << "Alpha: " << alpha << "; beta: " << beta << std::endl;
        ratio = (alpha + k)/(alpha + beta + n);
	//Rcpp::Rcout << "ratio: " << ratio << std::endl;
        rand1 = R::runif(0,1);
	// With probability 'ratio', add a 1. ie. if many 1s already, less likely
	// to add more 1s depending on alpha and beta
        if(rand1 < ratio){
          x(locs1(j)) = 1;
          k++;
        } else {
          x(locs1(j)) = 0;
        }
        n++;
      }
    } else {
      // Otherwise, swap the contents of N random locations
      maxI_indiv = x.size();
      IntegerVector moves;
      for(int j = 0; j < nInf; j++){
        id1 = floor(R::runif(0,1)*x.size());
        moveMax = moveSizes[indiv];
        move = floor(R::runif(0,1)*2*moveMax) - moveMax;
        id2 = id1 + move;
	while(id2 < 0) id2 += maxI_indiv;
        //if(id2 < 0) id2 = (move + id1) % maxI_indiv;
	if(id2 >= maxI_indiv) id2 = (id2 % maxI_indiv);
	tmp = x[id1];
        x[id1] = x[id2];
        x[id2] = tmp;
      }
    }
    newInfHist.submat(indiv, ageMask[indiv]-1, indiv, maxI-1) = x;
  }
  return(newInfHist);
}

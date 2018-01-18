#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

// [[Rcpp::export]]
arma::mat inf_hist_prop_cpp(arma::mat infHist, IntegerVector sampledIndivs, IntegerVector ageMask,IntegerVector moveSizes, IntegerVector nInfs,
			    double alpha, double beta, NumericVector randNs) {
  
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
  
  for(int i = 0; i < sampledIndivs.size(); ++i){
    indiv = sampledIndivs[i]-1;
    nInf = nInfs[indiv];
    x = newInfHist.submat(indiv, ageMask[indiv]-1, indiv, maxI-1);
    samps = seq_len(x.n_cols);
    if(randNs[i] < 1.0/2.0){
      locs = RcppArmadillo::sample(samps, nInf, FALSE, NumericVector::create());
      locs1 = as<arma::uvec>(locs)-1;
      y = newInfHist.row(indiv);
      y = y.elem(locs1);
      k = accu(x) - accu(y);
      n = x.size() - nInf;
  
      for(int j = 0; j < nInf; ++j){
        ratio = (alpha + k)/(alpha + beta + n);
        rand1 = R::runif(0,1);
        if(rand1 < ratio){
          x(locs1(j)) = 1;
          k++;
        } else {
          x(locs1(j)) = 0;
        }
        n++;
      }
    } else {
      maxI_indiv = x.size()-1;
      IntegerVector moves;
      for(int j = 0; j < nInf; j++){
        id1 = floor(R::runif(0,1)*x.size());
        moveMax = moveSizes[indiv];
        move = floor(R::runif(0,1)*2*moveMax) - moveMax;
        id2 = id1 + move;
        if(id2 < 0) id2 = maxI_indiv + id2;
        if(id2 > maxI_indiv) id2 = id2 - maxI_indiv;
        tmp = x[id1];
        x[id1] = x[id2];
        x[id2] = tmp;
      }
    }
    newInfHist.submat(indiv, ageMask[indiv]-1, indiv, maxI-1) = x;
  }
  return(newInfHist);
}

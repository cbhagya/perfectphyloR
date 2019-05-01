#include <RcppArmadillo.h>

//' Compute the dCor coefficient using Rcpp
//' 
//' @param mDx A numeric matrix of pairwise distances.
//' @param mDy A second numeric matrix of pairwise distances.
//' @param mC See the equation 2.4 in Josse and Homes manual.
//' 
//' @keywords Internal
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat dCor(arma::mat mDx, arma::mat mDy, arma::mat mC){
  
  arma::mat X = mC*mDx*mC;
  arma::mat Y = mC*mDy*mC;
  
  // Get trace(X'Y)
  arma::mat num1 = arma::diagvec(X*Y);
  // trace(X'Y) = trace(XY), since X is symmetric matrix.
  arma::mat num = sum(num1);
  
  // Compute ||X||*||Y||
  
  // 1. Compute ||X||
  // ||X|| = sqrt(trace(X*X')) = sqrt(trace(X*X))
  arma::mat deno1 = arma::diagvec(X*X);
  arma::mat deno1_1 = arma::sum(deno1);
  
  // 2. Compute ||Y||
  // ||Y|| = sqrt(trace(Y*Y')) = sqrt(trace(Y*Y))
  arma::mat deno2 = arma::diagvec(Y*Y);
  arma::mat deno2_1 = arma::sum(deno2);
  
  // Compute ||X||*||Y|| = sqrt(trace(XX')*trace(YY'))
  arma::mat deno = arma::sqrt(deno1_1*deno2_1);
  
  arma::mat dcor = num/deno ;
  // Equation 3.3 in Julie Josse, DOI: 10.1214/16-SS116
  return dcor;
  
 
}
/*** R
sourceCpp("dCorCoeff.cpp")
*/

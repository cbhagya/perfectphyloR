#include <RcppArmadillo.h>

//' Compute the RV coefficient using Rcpp
//' 
//' @param mDx A numeric matrix of pairwise distances.
//' @param mDy A second numeric matrix of pairwise distances.
//' @param mC See the equation 2.4 in Josse and Homes manual.
//' 
//' @keywords Internal
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat RVcoeff(arma::mat mDx, arma::mat mDy, arma::mat mC){
  
  arma::mat X = mC*arma::square(mDx)*mC;
  arma::mat Y = mC*arma::square(mDy)*mC;
  
  arma::mat num1 = arma::diagvec(X*Y);
  // tr(X'Y); X' = X since X is a symmetric matrix.
  arma::mat num = sum(num1); 
  // 
  arma::mat deno1 = arma::diagvec(X*X);
  // ||X|| = sqrt(tr(XX'))
  // Get tr(XX')
  arma::mat deno1_1 = arma::sum(deno1); 
  
  
  arma::mat deno2 = arma::diagvec(Y*Y);
  // ||Y|| = sqrt(tr(YY'))
  // Get tr(YY')
  arma::mat deno2_1 = arma::sum(deno2);
  
  // ||X||*||Y|| = sqrt(tr(XX') * tr(YY'))
  arma::mat deno = arma::sqrt(deno1_1*deno2_1);
  
  arma::mat RV = num/deno ;
  //double num = arma::as_scalar(nu);
  
  return RV;
  
  
}

/*** R
sourceCpp("RVcoeff.cpp")
 */

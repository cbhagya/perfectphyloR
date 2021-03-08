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

double dCor(arma::mat mDx, arma::mat mDy, arma::mat mC){
  
  arma::mat Wx = mC*mDx*mC; 
  arma::mat Wy = mC*mDy*mC;
  
  
  double num = arma::trace(Wx*Wy);
  
  // ||Wx|| = sqrt(trace(Wx*Wx))
  double deno1 = arma::trace(Wx*Wx);
  double deno2 = arma::trace(Wy*Wy);
  
  
  double deno = sqrt(deno1*deno2); 
  
  double dcor = sqrt(num/deno);
  return dcor;
  
 
}
/*** R
sourceCpp("dCorCoeff.cpp")
*/

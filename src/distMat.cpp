#include <RcppArmadillo.h>

//' 
//' This function computes the pair wise distances of tips according to the branching order using Rcpp.
//' 
//' @param len Length of tip set 1.
//' @param len2 Length of tip set 2.
//' @param nHaplos Number of haplotypes.
//' @param mD Distance matrix
//' 
//' @keywords internal
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat distMatRcpp(int len, arma::vec len2, int nHaplos, arma::mat mD){  
  
  int lenSplitTips = len;
  arma::vec lenSplitTips2 = len2;
  int nHaps = nHaplos;
  arma::mat distMat = mD;
  
  arma::mat matri = arma::zeros(nHaps, nHaps);
  
  int hapI = 0;
  int hapJ = 0;
  
  for(int i = 0; i < lenSplitTips; i++){
    for(int ii = 0; ii < lenSplitTips2[i]; ii++){
      
      hapI = hapI + 1; 
      
      for(int j = 0; j < lenSplitTips; j++){
        for(int jj = 0; jj < lenSplitTips2[j]; jj++){
          hapJ = hapJ+1;
          
          matri(hapI-1, hapJ-1) = distMat(i, j);
          
        }
      }
      hapJ = 0;
    }
  }
  return matri;
}


/*** R
sourceCpp("distMat.cpp")
  */

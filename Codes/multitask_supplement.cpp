#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// arma::mat get_B0(const arma::mat& X0tX0,
//                  const arma::mat& Ip,
//                  const arma::mat& X0tY0,
//                  const arma::mat& Gamma00,
//                  const arma::mat& Gamma01,
//                  const arma::mat& U0,
//                  const arma::mat& U1,
//                  const double N,
//                  const double K,
//                  const double rho00,
//                  const double rho01) { // only for one source!
//   
//   arma::mat lhs = (1.0 / (N*K)) * X0tX0 + (rho00+rho01) * Ip;
//   arma::mat rhs = (1.0 / (N*K)) * X0tY0 + rho00 * (Gamma00 - U0) + rho01 * (Gamma01 - U1);
//   
//   arma::mat B0 = arma::solve(lhs, rhs);
//   return B0;
// }

// // [[Rcpp::export]]
// arma::mat get_B1(const arma::mat& X1tX1,
//                  const arma::mat& Ip,
//                  const arma::mat& X1tY1,
//                  const arma::mat& Gamma1,
//                  const arma::mat& V1,
//                  const double N,
//                  const double K,
//                  const double rho1) {
//                      
//   arma::mat lhs = (1.0 / (N*K)) * X1tX1 + rho1 * Ip;
//   arma::mat rhs = (1.0 / (N*K)) * X1tY1 + rho1 * (Gamma1 - V1);
//   
//   arma::mat B1 = arma::solve(lhs, rhs);
//   return B1;
// }

// [[Rcpp::export]]
arma::vec update_row_block_cpp(const arma::vec& Aj, double Xj_norm2, double lambda) {
  double norm_Aj = arma::norm(Aj, 2);  // Euclidean norm (L2)
  
  if (norm_Aj > lambda) {
    return ((1.0 - lambda / norm_Aj) * Aj / Xj_norm2);
  } else {
    return arma::zeros<arma::vec>(Aj.n_elem);
  }
}

// [[Rcpp::export]]
arma::mat update_Gamma00(const arma::mat& B0,
                         const arma::mat& U0,
                         double lambda0,
                         double rho00,
                         const int p,
                         const int K) {
  arma::mat Z = B0 + U0;
  arma::mat Gamma00 = arma::zeros<arma::mat>(p, K);
  
  double lambda_ratio = lambda0 / rho00;
  
  for (arma::uword j = 0; j < p; ++j) {
    arma::vec Aj = Z.row(j).t();  // extract row as column vector
    Gamma00.row(j) = update_row_block_cpp(Aj, 1.0, lambda_ratio).t(); // transpose back to row
  }
  return Gamma00;
}

// [[Rcpp::export]]
arma::mat update_Gamma01(const arma::mat& B0,
                         const arma::mat& U1,
                         const arma::mat& Gamma1,
                         double lambda1,
                         double rho01,
                         const int p,
                         const int K){
  arma::mat Z = B0 + U1;
  arma::mat Gamma01 = arma::zeros<arma::mat>(p, K);
  double lambda_ratio = lambda1 / rho01;
  
  for (arma::uword j = 0; j < p; ++j) {
    arma::vec Aj = Z.row(j).t() - Gamma1.row(j).t();  // Aj = (B0 + U1 - Gamma1)[j, ]
    arma::vec update = update_row_block_cpp(Aj, 1.0, lambda_ratio);
    Gamma01.row(j) = Gamma1.row(j) + update.t();  // add to Gamma1 row
  }
  return Gamma01;
}

// [[Rcpp::export]]
arma::mat update_Gamma1(const arma::mat& B1,
                         const arma::mat& V1,
                         const arma::mat& Gamma01,
                         double lambda1,
                         double rho1,
                         const int p,
                         const int K){
  arma::mat Z = B1 + V1;
  arma::mat Gamma1 = arma::zeros<arma::mat>(p, K);
  double lambda_ratio = lambda1 / rho1;
  
  for (arma::uword j = 0; j < p; ++j) {
    arma::vec Aj = Z.row(j).t() - Gamma01.row(j).t();  // Aj = (B0 + U1 - Gamma1)[j, ]
    arma::vec update = update_row_block_cpp(Aj, 1.0, lambda_ratio);
    Gamma1.row(j) = Gamma01.row(j) + update.t();  // add to Gamma1 row
  }
  return Gamma1;
}
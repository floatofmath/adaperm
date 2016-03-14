// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix all_combinations_cpp(int n, int k){
  int nout, j, l;
  std::vector<bool> v(n);
  std::fill(v.begin(), v.end() - n + k, true);
  nout = Rf_choose(n,k);
  j = 0;
  l = 0;
  IntegerMatrix cs(k,nout);
   do {
     for (int i = 0; i < n; ++i) {
       if (v[i]) {
         cs(l,j) = i+1;
         l = l+1;
       }
     }
     j = j+1;
     l = 0;
   } while (std::prev_permutation(v.begin(), v.end()));
   return cs;
}

// [[Rcpp::export]]
LogicalMatrix all_reassignments_cpp(int n, int k){
  int nout, j;
  LogicalVector v(n);
  std::fill(v.begin(), v.end() - n + k, true);
  nout = Rf_choose(n,k);
  j = 0;
  LogicalMatrix cs(n,nout);
   do {
     cs(_,j) = v;
     j = j+1;
   } while (std::prev_permutation(v.begin(), v.end()));
   return cs;
}


// C++ program for space optimized Dynamic Programming
// Solution of Binomial Coefficient
// borrowed from http://www.geeksforgeeks.org/space-and-time-efficient-binomial-coefficient/
// Returns value of Binomial Coefficient C(n, k)
// int binomialCoeff(int n, int k)
// {
//     int res = 1;
 
//     // Since C(n, k) = C(n, n-k)
//     if ( k > n - k )
//         k = n - k;
 
//     // Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
//     for (int i = 0; i < k; ++i)
//     {
//         res *= (n - i);
//         res /= (i + 1);
//     }
 
//     return res;
// }


// mat combinations_rec(const int n, const int k, vec vs){
//   if(k == 1){
//     mat rmx(1,n);
//     rmx.row(0) = vs.t();
//     return rmx;
//   } else if(k == n) {
//     mat rmx(k,1);
//     rmx.col(0) = vs;
//     return rmx;
//   } else {
//     int nr = binomialCoeff(n-1,k-1);
//     mat v0mx(1,nr);
//     vec vts = vs.tail(n-1);
//     v0mx.fill(vs(0));
//     return( join_rows(join_cols(v0mx,combinations_rec(n - 1,k - 1,vts)),
//                       combinations_rec(n-1,k,vts)));
//   }
// }


// // [[Rcpp::export]]
// IntegerMatrix combinations_cpp(int n, int k){
//   vec vs = linspace<vec>(1, n, n);
//   mat cs = combinations_rec(n,k,vs);
//   return(wrap(cs));
// }

// // [[Rcpp::export]]
// NumericMatrix subsamples_cpp(NumericVector xs, const int k){
//   vec vs = as<vec>(xs);
//   int n = xs.size();
//   mat cs = combinations_rec(n,k,xs);
//   return(wrap(cs));
// }

// mat bincombinations_rec(const int k){
//   if(k == 1){
//     mat vmx =linspace<mat>(0,1,2);
//     return vmx.t();
//   } else {
//     int k2 = pow(2,k-1);
//     mat v0(1,k2);
//     mat v1(1,k2);
//     v0.fill(0);
//     v1.fill(1);
//     return( join_rows(join_cols(v0,bincombinations_rec(k-1)),
//                       join_cols(v1,bincombinations_rec(k-1))) );
//   }
// }

// // [[Rcpp::export]]
// IntegerMatrix bincombinations_cpp(const int p){
//   IntegerMatrix retval(pow(2,p),p);
//   int n;
//   for (n=1;n<p+1;n++) {
//     int p2 = pow(2,p);
//     int k2 = p2/pow(2,n);
//     IntegerVector v(2*k2);
//     std::fill(v.begin(),v.end(),1);
//     std::fill_n(v.begin(),k2,0);
//     retval(_, n-1) = rep(v,p2);
//   }
//   return retval;
// }


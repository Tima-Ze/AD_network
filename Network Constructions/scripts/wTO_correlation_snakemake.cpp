//#include <Rcpp.h>
#include <iostream>
#include <string>
#include <RcppArmadillo.h>
#include <algorithm> 

using namespace arma;
using namespace Rcpp;
using namespace std;


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

struct asset_info {
	double sum, sum2, stdev;
};

//[correlation matrix](http://en.wikipedia.org/wiki/Correlation_and_dependence).
// n,sX,sY,sXY,sX2,sY2
// cor = ( n * sXY - sX * sY ) / ( sqrt(n * sX2 - sX^2) * sqrt(n * sY2 - sY^2) )
inline asset_info compute_asset_info(const NumericMatrix& mat, 
	const int icol, const int rstart, const int rend) {
	double sum, sum2;
    sum = sum2 = 0;

	for (int r = rstart; r < rend; r++) {
		double d = mat(r, icol);
		sum += d;
		sum2 += pow(d,2);
	}
	
	asset_info res;
		res.sum = sum;
		res.sum2 = sum2;
		res.stdev = sqrt((rend-rstart) * sum2 - pow(sum, 2));
					
	return res;
}



/*------ Parallel version ------*/

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

// pre-compute sum and stdev
struct cor_p1 : public Worker {
	const RMatrix<double> mat;
	const int rstart, rend, nperiod;

	RVector<double> rsum, rstdev;

	cor_p1(const NumericMatrix& mat, const int rstart, const int rend,
		NumericVector rsum, NumericVector rstdev)
		: mat(mat), rstart(rstart), rend(rend), nperiod(rend - rstart), rsum(rsum), rstdev(rstdev) { }

	void operator()(size_t begin, size_t end) {
		for (size_t c = begin; c < end; c++) {
			double sum, sum2;
			sum = sum2 = 0;

			for (int r = rstart; r < rend; r++) {
				double d = mat(r,c);
				sum += d;
				sum2 += pow(d,2);
			}

			rsum[c] = sum;
			rstdev[c] = sqrt(nperiod * sum2 - pow(sum,2));
		}
	}
};
// compute correlation
struct cor_p2 : public Worker {
	const RMatrix<double> mat;
	const int rstart, rend, nperiod;
	const RVector<double> sum, stdev;
    
	 RMatrix<double> rmat;

	cor_p2(const NumericMatrix& mat, const int rstart, const int rend,
		const NumericVector& sum, const NumericVector& stdev, 
		NumericMatrix rmat)
		: mat(mat), rstart(rstart), rend(rend), nperiod(rend - rstart), sum(sum), stdev(stdev), rmat(rmat) {}

	void operator()(size_t begin, size_t end) {
		for (size_t c1 = begin; c1 < end; c1++) {
			for (size_t c2 = 0; c2 < c1; c2++) {
				double sXY = 0;
				for (int r = rstart; r < rend; r++)
					sXY += mat(r,c1) * mat(r,c2);

				rmat(c1,c2) = (nperiod * sXY - sum[c1] * sum[c2]) / (stdev[c1] * stdev[c2]);         
			}
		}
	}
};

inline NumericMatrix cp_cor_helper(const NumericMatrix& mat, const int rstart, const int rend) {
	int nc = mat.ncol();
	NumericVector rsum(nc), rstdev(nc);

	cor_p1 p1(mat, rstart, rend, rsum, rstdev);
	parallelFor(0, nc, p1);
	
	NumericMatrix rmat(nc, nc);
	
	cor_p2 p2(mat, rstart, rend, rsum, rstdev, rmat);
	parallelFor(0, nc, p2);
	//std::cout << rmat << std::endl;
	return rmat;
}

inline NumericMatrix cp_cor1(NumericMatrix vmat) {
	//std::cout << vmat << std::endl;
	int ll = vmat.ncol();
	int pp = vmat.nrow();
	for(int i = 0; i < ll; ++i){
		for(int j = 0; j < i; ++j){
			vmat[i*pp+j] = vmat[j*pp+i];
		}
	}    
	return vmat;
}

// [[Rcpp::export]]
NumericMatrix cp_cor(NumericMatrix mat) {
	//std::cout << rmat << std::endl;
	NumericMatrix rmat = cp_cor_helper(mat, 0, mat.nrow());
	//std::cout << rmat << std::endl;
	return cp_cor1(rmat);
	//return rmat;
}


// [[Rcpp::export]]
arma::uvec calc_ranks(const arma::vec& da_ta) {
  return (arma::sort_index(arma::sort_index(da_ta)) + 1);
}  // end calc_ranks

arma::uvec calc_ranks1(const arma::vec& iwas) {
  return (arma::sort_index(arma::sort_index(iwas)) + 1);
}  // end calc_ranks


// [[Rcpp::export]]
NumericMatrix RankMatrix(NumericMatrix cmat) {
	int g = cmat.ncol();

	for(int i = 0; i < g; ++i){
	  NumericVector p = cmat.column(i);
	  p =  calc_ranks1(p);
	  cmat.column(i) = p;
		}
	return cmat;
}

// [[Rcpp::export]]
NumericMatrix spearman_cor(NumericMatrix mat) {
  NumericMatrix smat = RankMatrix(mat);
  smat = cp_cor_helper(smat, 0, mat.nrow());
  //std::cout << rmat << std::endl;
  return cp_cor1(smat);
}





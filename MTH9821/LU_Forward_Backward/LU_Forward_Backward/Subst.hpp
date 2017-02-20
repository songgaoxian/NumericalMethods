//define functions forward substitution and backward substitution

#ifndef Subst_HPP
#define Subst_HPP

#include"MatrixLU.hpp"
#include"BandMatrixLU.hpp"

//forward_substitution
using Vector = Eigen::VectorXd;
//for regular matrix
Vector forward_subst(Matrix L, Vector b) {
	Vector result(L.cols());
	result(0) = b(0) / L(0, 0);
	double sum = 0;
	for (int k = 1; k < L.rows(); ++k) {
		sum = 0;
		for (int i = 0; i < k; ++i) {
			sum += L(k, i)*result(i);
		}
		result(k) = (b(k) - sum) / L(k, k);
	}
	return result;
}
//for band matrix
Vector forward_subst_bandm(BandMatrix L, Vector b) {
	int m = L.subs();//get bandwidth
	Vector result(L.cols());
	result(0) = b(0) / item(L,0, 0);
	double sum = 0;
	for (int k = 1; k < L.rows(); ++k) {
		sum = 0;
		for (int i = std::max(0, k - m); i < k; ++i) {
			sum += item(L,k, i)*result(i);
		}
		result(k) = (b(k) - sum) / item(L,k, k);
	}
	return result;
}
//bandwith=1
Vector forward_subst_band1(BandMatrix L, Vector b) {
	Vector result(L.cols());
	result(0) = b(0) / item(L,0, 0);
	double sum;
	for (int k = 1; k < L.rows(); ++k) {
		sum = item(L,k, k-1)*result(k-1);
		result(k) = (b(k) - sum) / item(L,k, k);
	}
	return result;
}

Vector backward_subst(Matrix U, Vector b) {
	Vector result(U.cols());
	int c = U.cols();
	result(c - 1) = b(c - 1) / U(c - 1, c - 1);
	double sum = 0;
	for (int k = c - 2; k >= 0; --k) {
		sum = 0;
		for (int i = c - 1; i > k; --i) {
			sum += result(i)*U(k, i);
		}
		result(k) = (b(k) - sum) / U(k, k);
	}
	return result;
}
//for bandmatrix
Vector backward_subst_bandm(BandMatrix U, Vector b) {
	int m = U.supers();
	Vector result(U.cols());
	int c = U.cols();
	result(c - 1) = b(c - 1) / item(U,c - 1, c - 1);
	double sum = 0;
	for (int k = c - 2; k >= 0; --k) {
		sum = 0;
		for (int i = std::min(c - 1,k+m); i > k; --i) {
			sum += result(i)*item(U,k, i);
		}
		result(k) = (b(k) - sum) / item(U,k, k);
	}
	return result;
}
//bandwidth=1
Vector backward_subst_band1(BandMatrix U, Vector b) {
	Vector result(U.cols());
	int c = U.cols();
	result(c - 1) = b(c - 1) / item(U,c - 1, c - 1);
	double sum;
	for (int k = c - 2; k >= 0; --k) {
		sum = result(k+1)*item(U,k, k+1);
		result(k) = (b(k) - sum) / item(U,k, k);
	}
	return result;
}

Vector diagonal_subst(BandMatrix D, Vector b) {
	Vector result(D.cols());
	for (int i = 0; i < D.cols(); ++i) {
		result(i) = b(i) / item(D,i, i);
	}
	return result;
}

Vector diagonal_subst(Matrix D, Vector b) {
	Vector result(D.cols());
	for (int i = 0; i < D.cols(); ++i) {
		result(i) = b(i) / D(i, i);
	}
	return result;
}

#endif
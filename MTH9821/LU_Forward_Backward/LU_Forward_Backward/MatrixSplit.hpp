//split matrixes
#ifndef MatrixSplit_HPP
#define MatrixSplit_HPP

#include"BandMatrixLU.hpp"
#include"MatrixLU.hpp"
#include"Subst.hpp"

std::tuple<Matrix, Matrix, Matrix> SplitMatrix(Matrix A) {
	Matrix D(A.rows(),A.cols());
	D.setZero();
	for (int i = 0; i < A.cols(); ++i) {
		D(i, i) = A(i, i);
	}
	Matrix L(A.rows(), A.cols());
	L.setZero();
	for (int i = 0; i < A.rows(); ++i) {
		for (int j = 0; j < i; ++j) {
			L(i, j) = A(i, j);
		}
	}
	Matrix U(A.rows(), A.cols());
	U.setZero();
	U = A - D - L;
	std::tuple<Matrix, Matrix, Matrix> result(D, L, U);
	return result;
}
std::tuple<BandMatrix, BandMatrix, BandMatrix> SplitMatrix(BandMatrix A) {
	BandMatrix D(A.rows(), A.cols(), 0, 0);
	BandMatrix L(A.rows(), A.cols(), 0, A.subs());
	BandMatrix U(A.rows(), A.cols(), A.supers(), 0);
	for (int i = 0; i < A.rows(); ++i) {
		item(D, i, i) = item(A, i, i);
		item(L, i, i) = 0;
		item(U, i, i) = 0;
	}
	
	for (int i = 0; i < A.cols(); ++i) {
		for (int j = i+1; j <= std::min(A.rows()-1,i+A.subs()); ++j) {
			item(L, j, i) = item(A, j, i );
		}
		for (int j = i-1; j >=std::max(0,i-A.supers()); --j) {
			item(U, j, i) = item(A, j, i);
		}
	}
	std::tuple<BandMatrix, BandMatrix, BandMatrix> result(D, L, U);
	return result;
}

double Norm2(Vector b) {
	double result = 0;
	for (int i = 0; i < b.rows(); ++i) {
		result += b(i)*b(i);
	}
	return std::sqrt(result);
}


#endif
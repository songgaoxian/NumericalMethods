//implement Cholesky decomposition

#ifndef Cholesky_HPP
#define Cholesky_HPP

#include"BandMatrixLU.hpp"
#include"MatrixLU.hpp"
//regular matrix
Matrix Cholesky_matrix(Matrix m) {
	Matrix U(m.rows(), m.cols());
	U.setZero();//initialize U as zero matrix
	//try Cholesky decomposition
	try{
		for (int i = 0; i < m.rows()-1; ++i) {
			U(i, i) = std::sqrt(m(i, i));
			//get values for ith row
			for (int j = i + 1; j < m.cols(); ++j) {
				U(i, j) = m(i, j) / U(i, i);
			}
			//update m
			for(int j=i+1;j<m.rows();++j)
				for (int k = j; k < m.cols(); ++k) {
					m(j, k) = m(j, k) - U(i, j)*U(i, k);
					m(k, j) = m(j, k);
				}
		}
		int size = m.rows();
		U(size - 1, size - 1) = std::sqrt(m(size - 1, size - 1));
		return U;
	}
	catch (...) {
		std::cout << "Given matrix is not spd\n";
	}
}
//bandmatrix 
BandMatrix Cholesky_band(BandMatrix m) {
	int band = m.subs();
	BandMatrix U(m.rows(), m.cols(),band,0);
	//try Cholesky decomposition
	try {
		for (int i = 0; i < m.rows() - 1; ++i) {
			item(U,i, i) = std::sqrt(item(m,i, i));
			//get values for ith row
			for (int j = i + 1; j < std::min(m.cols(),i+band+1); ++j) {
				item(U,i, j) = item(m,i, j) / item(U,i, i);
			}
			//update m
			for (int j = i + 1; j<std::min(m.rows(),i+band+1); ++j)
				for (int k = j; k < std::min(m.cols(),i+band+1); ++k) {
					item(m,j, k) = item(m,j, k) - item(U,i, j)*item(U,i, k);
					item(m,k, j) = item(m,j, k);
				}
		}
		int size = m.rows();
		item(U,size - 1, size - 1) = std::sqrt(item(m,size - 1, size - 1));
		return U;
	}
	catch (...) {
		std::cout << "Given matrix is not spd\n";
	}
}
//tridiagonal
BandMatrix Cholesky_band1(BandMatrix m) {
	BandMatrix U(m.rows(), m.cols(), 1, 0);
	//try Cholesky decomposition
	try {
		for (int i = 0; i < m.rows() - 1; ++i) {
			item(U, i, i) = std::sqrt(item(m, i, i));
			//get values for ith row
			item(U, i, i+1) = item(m, i, i+1) / item(U, i, i);
			//update m
			item(m, i+1, i+1) = item(m, i+1, i+1) - item(U, i, i+1)*item(U, i, i+1);
		}
		int size = m.rows();
		item(U, size - 1, size - 1) = std::sqrt(item(m, size - 1, size - 1));
		return U;
	}
	catch (...) {
		std::cout << "Given matrix is not spd\n";
	}
}

#endif
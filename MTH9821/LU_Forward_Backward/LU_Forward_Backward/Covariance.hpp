//compute covariance
#ifndef Covariance_HPP
#define Covariance_HPP
#include"LinearSolver.hpp"

Matrix GetCov(Matrix X) {
	int row = X.rows(), col = X.cols();
	Matrix Cov(col, col);
	Vector Mean(col);
	for (int i = 0; i < col; ++i) {
		double sum = 0.0;
		for (int j = 0; j < row; ++j) {
			sum += X(j, i);
		}
		Mean(i) = sum / double(row);
	}
	for (int prod1 = 0; prod1 < col; ++prod1) {
		for (int prod2 = prod1; prod2 < col; ++prod2) {
			double sum = 0.0;
			for (int i = 0; i < row; ++i) {
				sum += (X(i, prod1) - Mean(prod1))*(X(i, prod2) - Mean(prod2));
			}
			Cov(prod1, prod2) = sum / double(row - 1.0);
			Cov(prod2, prod1) = Cov(prod1, prod2);
		}
	}
	return Cov;
}

Vector LinearRegressionCoef(Vector y, Matrix X) {
	return lu_no_pivoting_solver(X.transpose()*X, X.transpose()*y);
}


#endif
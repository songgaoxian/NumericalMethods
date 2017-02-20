//implement linear solver for Cholesky Decomposition and LU Decomposition
#ifndef LinearSolver_HPP
#define LinearSolver_HPP

#include"CholeskyDecomposition.hpp"
#include"Subst.hpp"

Vector Linear_chol_solver(Matrix A, Vector b) {
	Vector x(A.cols());//store result
	Matrix U(A.rows(), A.cols()),Ut(A.rows(),A.cols());
	U = Cholesky_matrix(A);
	Vector y(A.cols());//intermediate result
	Ut = U.transpose();
	y = forward_subst(Ut, b);
	x = backward_subst(U, y);
	//std::cout << U << '\n';
	return x;
}

Vector Linear_chol_solver(BandMatrix A, Vector b) {
	Vector x(A.cols()), y(A.cols());
	BandMatrix U(A.rows(), A.cols(), A.supers(), 0), Ut(A.rows(), A.cols(), 0, A.subs());
	U = Cholesky_band(A);
	Ut = Transpose(U);
	y = forward_subst_bandm(Ut, b);
	x = backward_subst_bandm(U, y);
	//std::cout << U.toDenseMatrix() << '\n';
	return x;
}
//tridiagonal
Vector Linear_tri_solver(BandMatrix A, Vector b) {
	Vector x(A.cols()), y(A.cols());
	BandMatrix U(A.rows(), A.cols(), A.supers(), 0), Ut(A.rows(), A.cols(), 0, A.subs());
	U = Cholesky_band1(A);
	Ut = Transpose(U);
	y = forward_subst_band1(Ut, b);
	x = backward_subst_band1(U, y);
	//std::cout << U.toDenseMatrix() << '\n';
	return x;
}

#endif
//use lu decomposition to solve linear system and inverse computation
#ifndef LU_Solver_HPP
#include"MatrixLU.hpp"
#include"Subst.hpp"

Vector lu_no_pivoting_solver(Matrix A, Vector b) {
	std::tuple<Matrix, Matrix> result = lu_no_pivoting(A);
	Matrix L(std::get<0>(result)), U(std::get<1>(result));
	Vector y(forward_subst(L, b));
	Vector x(backward_subst(U, y));
	return x;
}
Vector lu_no_pivoting_solver(BandMatrix A, Vector b) {
	std::tuple<BandMatrix, BandMatrix> result = lu_no_pivoting(A);
	BandMatrix L(std::get<0>(result)), U(std::get<1>(result));
	Vector y(forward_subst_bandm(L, b));
	Vector x(backward_subst_bandm(U, y));
	return x;
}
Vector lu_row_pivoting_solver(Matrix A, Vector b) {
	std::tuple<Matrix, Matrix,Matrix> result = lu_row_pivoting(A);
	Matrix P(std::get<0>(result)), L(std::get<1>(result)), U(std::get<2>(result));
	b = P*b;
	Vector y(forward_subst(L, b));
	Vector x(backward_subst(U, y));
	return x;
}
Matrix inverse_regular(Matrix A) {
	Matrix I(A.rows(), A.cols());
	I.setIdentity();
	Matrix Result(A.rows(), A.cols());
	for (int i = 0; i < A.rows(); ++i) {
		Result.col(i) = lu_row_pivoting_solver(A, I.col(i));
	}
	return Result;
}
Matrix inverse_triangle(Matrix A,bool isUpper) {
	Matrix I(A.rows(), A.cols());
	I.setIdentity();
	Matrix Result(A.rows(), A.cols());
	if (isUpper) {
		for (int i = 0; i < A.rows(); ++i) {
			Result.col(i) = backward_subst(A, I.col(i));
		}
	}
	else {
		for (int i = 0; i < A.rows(); ++i) {
			Result.col(i) = forward_subst(A, I.col(i));
		}
	}
	return Result;
}



#endif
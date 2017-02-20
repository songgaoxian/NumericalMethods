//implement iteration method
#ifndef Iteration_HPP
#define Iteration_HPP

#include"MatrixLU.hpp"
#include"BandMatrixLU.hpp"
#include"Subst.hpp"
#include"MatrixSplit.hpp"
#include<fstream>
#include<iomanip>

void WriteToCsv(std::string filename, Matrix A) {
	std::ofstream file(filename);
	for (int i = 0; i < A.rows(); ++i) {
		for (int j = 0; j < A.cols(); ++j) {
			file << std::setprecision(18) << A(i, j) << ",";
		}
		file << "\n";
	}
}
void WriteToCsv(std::string filename, Vector b, double r) {
	std::ofstream file(filename);
	for (int i = 0; i < b.rows(); ++i) {
		file << std::setprecision(18) << b(i) << "\n";
	}
	file << "residualORcount\n";
	file << std::setprecision(18) << r << "\n";
}

Vector vectorMultBand(BandMatrix A, Vector X) {
	Vector result(A.rows());		
		result.setZero();
		for (int i = 0; i < A.rows(); ++i) {
			for (int j = std::max(0, i - A.subs()); j <= std::min(A.cols() - 1, i + A.supers()); ++j) {
				result(i) += item(A, i, j)*X(j);
			}
		}
	return result;
}
BandMatrix scalarMultBand(BandMatrix A, double w) {
		for (int i = 0; i < A.rows(); ++i) {
			for (int j = std::max(0, i - A.subs()); j <= std::min(A.cols() - 1, i + A.supers()); ++j) {
				item(A, i, j) *= w;
			}
		}
	
	return A;
}
BandMatrix operator+(BandMatrix A, BandMatrix B) {
	BandMatrix result(A.rows(), A.cols(), std::max(A.supers(), B.supers()), std::max(A.subs(), B.subs()));
		int sub = std::min(A.subs(), B.subs()), super = std::min(A.supers(), B.supers());
		for (int i = 0; i < A.rows(); ++i) {
			for (int j = std::max(0, i - A.subs()); j <= std::min(A.cols() - 1, i + A.supers()); ++j) {
				item(result, i, j) = item(A, i, j);
			}
		}
		for (int i = 0; i < A.rows(); ++i) {
			for (int j = std::max(0, i - B.subs()); j <= std::min(A.cols() - 1, i + B.supers()); ++j) {
				item(result, i, j) = item(B, i, j);
			}
		}
		for (int i = 0; i < A.rows(); ++i) {
			for (int j = std::max(0, i - sub); j <= std::min(A.cols() - 1, i + super); ++j) {
					item(result, i, j) = item(A, i, j) + item(B, i, j);
			}
		}
	return result;
}

std::tuple<Vector,int> Jacobi_residual(Matrix A, Vector b, double tol, Vector x0) {
	std::tuple<Matrix, Matrix, Matrix> split_result = SplitMatrix(A);
	Matrix D(A.rows(),A.cols()), L(A.rows(),A.cols()), U(A.rows(),A.cols());
	D = std::get<0>(split_result);
	L = std::get<1>(split_result);
	U = std::get<2>(split_result);
	Vector x(x0), r(A.rows());
	r = b - A*x;
	double stop_err = tol*Norm2(r);
	Vector b_new(A.rows());
	b_new = diagonal_subst(D, b);
	int ic = 0;
	Vector temp(A.rows());
	while (Norm2(r) > stop_err) {
		temp = L*x + U*x;
		x = -diagonal_subst(D, temp) + b_new;
		r = b - A*x;
		ic++;
		if (ic == 1) WriteToCsv("jaco1.csv", x);
		if (ic == 2) WriteToCsv("jaco2.csv", x);
		if (ic == 3) WriteToCsv("jaco3.csv", x);
	}
	std::tuple<Vector, int> result(x, ic);
	return result;
}
std::tuple<Vector, int> Jacobi_residual(BandMatrix A, Vector b, double tol,Vector x0) {
	std::tuple<BandMatrix, BandMatrix, BandMatrix> split_result = SplitMatrix(A);
	BandMatrix D(std::get<0>(split_result)), L(std::get<1>(split_result)), U(std::get<2>(split_result));
	Vector x(x0), r(A.rows());
	r = b - vectorMultBand(A,x);
	double stop_err = tol*Norm2(r);
	Vector b_new(A.rows());
	b_new = diagonal_subst(D, b);
	int ic = 0;
	Vector temp(A.rows());
	while (Norm2(r) > stop_err) {
		temp = vectorMultBand(L,x) + vectorMultBand(U,x);
		x = -diagonal_subst(D, temp) + b_new;
		r = b - vectorMultBand(A,x);
		ic++;
	}
	std::tuple<Vector, int> result(x, ic);
	return result;
}
std::tuple<Vector, int> Jacobi_consecutive(Matrix A, Vector b, double tol, Vector x0) {
	std::tuple<Matrix, Matrix, Matrix> split_result = SplitMatrix(A);
	Matrix D(A.rows(), A.cols()), L(A.rows(), A.cols()), U(A.rows(), A.cols());
	D = std::get<0>(split_result);
	L = std::get<1>(split_result);
	U = std::get<2>(split_result);
	Vector xold(x0), xnew(x0.rows()), diff(x0.rows());
	Vector b_new(A.rows());
	b_new = diagonal_subst(D, b);
	int ic = 0;
	Vector temp(A.rows());
	do {
		temp = L*xold + U*xold;
		xnew = -diagonal_subst(D, temp) + b_new;
		diff = xnew - xold;
		xold = xnew;
		ic++;
	} while (Norm2(diff) > tol);
	std::tuple<Vector, int> result(xnew, ic);
	return result;
}
std::tuple<Vector, int> Jacobi_consecutive(BandMatrix A, Vector b, double tol, Vector x0) {
	std::tuple<BandMatrix, BandMatrix, BandMatrix> split_result = SplitMatrix(A);
	BandMatrix D(std::get<0>(split_result)), L(std::get<1>(split_result)), U(std::get<2>(split_result));
	Vector xold(x0), xnew(x0.rows()), diff(x0.rows());
	Vector b_new(A.rows());
	b_new = diagonal_subst(D, b);
	int ic = 0;
	Vector temp(A.rows());
	do {
		temp = vectorMultBand(L, xold) + vectorMultBand(U, xold);
		xnew = -diagonal_subst(D, temp) + b_new;
		diff =xnew-xold;
		xold = xnew;
		ic++;
	} while (Norm2(diff) > tol);
	std::tuple<Vector, int> result(xnew, ic);
	return result;
}

std::tuple<Vector, int> GaussSiedel_residual(Matrix A, Vector b, double tol,Vector x0) {
	std::tuple<Matrix, Matrix, Matrix> split_result = SplitMatrix(A);
	Matrix D(std::get<0>(split_result)), L(std::get<1>(split_result)), U(std::get<2>(split_result));
	Vector x(x0), r(A.rows());
	r = b - A*x;
	double stop_err = tol*Norm2(r);
	Vector b_new(A.rows());
	Matrix DL(L + D);
	b_new = forward_subst(DL, b);
	int ic = 0;
	Vector temp(A.rows());
	while (Norm2(r) > stop_err) {
		temp = U*x;
		x = -forward_subst(DL, temp) + b_new;
		ic++;
		r = b - A*x;
		if (ic == 1) WriteToCsv("gauss1.csv", x);
		if (ic == 2) WriteToCsv("gauss2.csv", x);
		if (ic == 3) WriteToCsv("gauss3.csv", x);
	}
	std::tuple<Vector, int> result(x, ic);
	return result;
}

std::tuple<Vector, int> GaussSiedel_residual(BandMatrix A, Vector b, double tol,Vector x0) {
	std::tuple<BandMatrix, BandMatrix, BandMatrix> split_result = SplitMatrix(A);
	BandMatrix D(std::get<0>(split_result)), L(std::get<1>(split_result)), U(std::get<2>(split_result));
	Vector x(x0), r(A.rows());
	r = b - vectorMultBand(A,x);
	double stop_err = tol*Norm2(r);
	Vector b_new(A.rows());
	BandMatrix DL(L + D);
	b_new = forward_subst_bandm(DL, b);
	int ic = 0;
	Vector temp(A.rows());
	while (Norm2(r) > stop_err) {
		temp = vectorMultBand(U,x);
		x = -forward_subst_bandm(DL, temp) + b_new;
		ic++;
		r = b - vectorMultBand(A,x);
	}
	std::tuple<Vector, int> result(x, ic);
	return result;
}

//consecutive criteria
std::tuple<Vector, int> GaussSiedel_consecutive(Matrix A, Vector b, double tol, Vector x0) {
	std::tuple<Matrix, Matrix, Matrix> split_result = SplitMatrix(A);
	Matrix D(std::get<0>(split_result)), L(std::get<1>(split_result)), U(std::get<2>(split_result));
	Vector xold(x0), diff(x0.rows()), xnew(x0.rows());
	Vector b_new(A.rows());
	Matrix DL(L + D);
	b_new = forward_subst(DL, b);
	int ic = 0;
	Vector temp(A.rows());
	do {
		temp = U*xold;
		xnew = -forward_subst(DL, temp) + b_new;
		diff = xnew - xold;
		ic++;
		xold = xnew;
	} while (Norm2(diff) > tol);
	std::tuple<Vector, int> result(xnew, ic);
	return result;
}

std::tuple<Vector, int> GaussSiedel_consecutive(BandMatrix A, Vector b, double tol, Vector x0) {
	std::tuple<BandMatrix, BandMatrix, BandMatrix> split_result = SplitMatrix(A);
	BandMatrix D(std::get<0>(split_result)), L(std::get<1>(split_result)), U(std::get<2>(split_result));
	Vector xold(x0), diff(x0.rows()), xnew(x0.rows());
	Vector b_new(A.rows());
	BandMatrix DL(L + D);
	b_new = forward_subst_bandm(DL, b);
	int ic = 0;
	Vector temp(A.rows());
	do {
		temp = vectorMultBand(U, xold);
		xnew = -forward_subst_bandm(DL, temp) + b_new;
		diff = xnew - xold;
		ic++;
		xold = xnew;
	} while (Norm2(diff) > tol);
	std::tuple<Vector, int> result(xnew, ic);
	return result;
}

std::tuple<Vector, int> SOR_residual(Matrix A, Vector b, double tol, double w,Vector x0) {
	std::tuple<Matrix, Matrix, Matrix> split_result = SplitMatrix(A);
	Matrix D(std::get<0>(split_result)), L(std::get<1>(split_result)), U(std::get<2>(split_result));
	Vector x(x0), r(A.rows());
	r = b - A*x;
	double stop_err = tol*Norm2(r);
	Vector b_new(A.rows());
	Matrix DwL(D + w*L);
	b_new = w*forward_subst(DwL, b);
	int ic = 0;
	Vector temp(A.rows());
	while (Norm2(r) > stop_err) {
		temp = ((1 - w)*D - w*U)*x;
		x = forward_subst(DwL, temp) + b_new;
		r = b - A*x;
		ic++;
		//need to check the criterion
		if (w < 1.0) {
			if (ic == 1) WriteToCsv("sor11.csv", x);
			if (ic == 2) WriteToCsv("sor12.csv", x);
			if (ic == 3) WriteToCsv("sor13.csv", x);
		}
		else {
			if (ic == 1) WriteToCsv("sor21.csv", x);
			if (ic == 2) WriteToCsv("sor22.csv", x);
			if (ic == 3) WriteToCsv("sor23.csv", x);
		}
	}
	std::tuple<Vector, int> result(x, ic);
	return result;
}

std::tuple<Vector, int> SOR_residual(BandMatrix A, Vector b, double tol, double w,Vector x0) {
	std::tuple<BandMatrix, BandMatrix, BandMatrix> split_result = SplitMatrix(A);
	BandMatrix D(std::get<0>(split_result)), L(std::get<1>(split_result)), U(std::get<2>(split_result));
	Vector x(x0),r(A.rows());
	r = b - vectorMultBand(A,x);
	double stop_err = tol*Norm2(r);
	Vector b_new(A.rows());
	BandMatrix DwL(D + scalarMultBand(L,w));
	b_new = w*forward_subst_bandm(DwL, b);
	int ic = 0;
	Vector temp(A.rows());
	while (Norm2(r) > stop_err) 
	{
		BandMatrix temp1(scalarMultBand(D,1-w)+ scalarMultBand(U,-w));
		temp = vectorMultBand(temp1, x);
		x = forward_subst_bandm(DwL, temp) + b_new;
		r = b - vectorMultBand(A,x);
		ic++;
	} 
	std::tuple<Vector, int> result(x, ic);
	return result;
}

std::tuple<Vector, int> SOR_consecutive(Matrix A, Vector b, double tol, double w, Vector x0) {
	std::tuple<Matrix, Matrix, Matrix> split_result = SplitMatrix(A);
	Matrix D(std::get<0>(split_result)), L(std::get<1>(split_result)), U(std::get<2>(split_result));
	Vector xold(x0), xnew(x0.rows()), diff(x0.rows());
	Vector b_new(A.rows());
	Matrix DwL(D + w*L);
	b_new = w*forward_subst(DwL, b);
	int ic = 0;
	Vector temp(A.rows());
	do {
		temp = ((1 - w)*D - w*U)*xold;
		xnew = forward_subst(DwL, temp) + b_new;
		diff = xnew - xold;
		xold = xnew;
		ic++;
	} while (Norm2(diff) > tol);
	std::tuple<Vector, int> result(xnew, ic);
	return result;
}

std::tuple<Vector, int> SOR_consecutive(BandMatrix A, Vector b, double tol, double w, Vector x0) {
	std::tuple<BandMatrix, BandMatrix, BandMatrix> split_result = SplitMatrix(A);
	BandMatrix D(std::get<0>(split_result)), L(std::get<1>(split_result)), U(std::get<2>(split_result));
	Vector xold(x0), xnew(x0.rows()), diff(x0.rows());
	Vector b_new(A.rows());
	BandMatrix DwL(D + scalarMultBand(L, w));
	b_new = w*forward_subst_bandm(DwL, b);
	int ic = 0;
	Vector temp(A.rows());
	do {
		BandMatrix temp1(scalarMultBand(D, 1 - w) + scalarMultBand(U, -w));
		temp = vectorMultBand(temp1, xold);
		xnew = forward_subst_bandm(DwL, temp) + b_new;
		diff = xnew - xold;
		xold = xnew;
		ic++;
	} while (Norm2(diff) > tol);
	std::tuple<Vector, int> result(xnew, ic);
	return result;
}



#endif
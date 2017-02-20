#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<vector>
#include<iomanip>
#include<boost\algorithm\string.hpp>
#include"LinearSolver.hpp"
#include"Iteration.hpp"
#include"LU_Solver.hpp"
#include"Covariance.hpp"
#include"Cubic_Spline.hpp"

void IniMatrixCsv(std::ifstream& file, Matrix& A) {
	std::string line, field;
	std::vector<std::string> dataline;
	while (std::getline(file, line)) {
		dataline.push_back(line);
	}

	std::stringstream is;
	for (int i = 0; i < dataline.size(); ++i) {
		std::vector<std::string> data;
		boost::split(data, dataline[i], boost::is_any_of(" ,\t"));
		for (int j = 0; j < data.size(); ++j) {
			is = std::stringstream(data[j]);
			is >> A(i, j);
		}
	}
}
void IniBandMatrixCsv(std::ifstream& file, BandMatrix& A) {
	std::string line, field;
	std::vector<std::string> dataline;
	while (std::getline(file, line)) {
		dataline.push_back(line);
	}
	double x;
	std::stringstream is;
	for (int i = 0; i < dataline.size(); ++i) {
		std::vector<std::string> data;
		boost::split(data, dataline[i], boost::is_any_of(" ,\t"));
		for (int j = 0; j < data.size(); ++j) {
			is = std::stringstream(data[j]);
			is >>x;
			if (x != 0) { item(A,i, j) = x; }
		}
	}
}
void IniVectorCsv(std::ifstream& file, Vector& A) {
	std::string line, field;
	std::vector<std::string> dataline;
	std::stringstream is;
	int i = 0;
	while (std::getline(file, line)) {
		is = std::stringstream(line);
		is >> A(i);
		i++;
	}
}


int main() {
	int N = 8;
	int row = 5, col =2,sub=1,super=1;
	std::cout.precision(18);
	double tol = 1E-6;
	Matrix discount(row, col);
	std::ifstream file("data.csv");
	IniMatrixCsv(file, discount);
	Matrix Rate(DiscountToRate(discount, 0.007));
	WriteToCsv("rate.csv", Rate);
	Matrix theM(std::get<0>(GetTridiagonalSys(Rate)).toDenseMatrix());
	WriteToCsv("theM.csv", theM);
	std::ifstream file2("data1.csv");
	Matrix time_cf(4, 2);
	IniMatrixCsv(file2, time_cf);
	double p = Price(time_cf, Rate);
	WriteToCsv("price.csv",discount.col(0), p);
	/*Vector K(rawData.col(0));
	Vector C(0.5*rawData.col(1) + 0.5*rawData.col(2));
	Vector P(0.5*rawData.col(3) + 0.5*rawData.col(4));
	Vector y(C - P);
	Matrix A(row, 2);
	A.col(0).setOnes();
	A.col(1) = -1.0*K;
	Vector x(LinearRegressionCoef(y, A));
	WriteToCsv("A.csv", A);
	WriteToCsv("y.csv", y);
	WriteToCsv("pvf_disc.csv", x);*/
	//std::ifstream file("data.csv");
	//IniMatrixCsv(file, discount);
	//Matrix Rate(DiscountToRate(discount, 0.0075));
	//Matrix Time_CF(row - 1, col);
	//std::ifstream file2("data1.csv");
	//IniMatrixCsv(file2, Time_CF);
	//double p = Price(Time_CF, Rate);
	//std::cout << p << '\n';
	//A4.setZero();
	/*for (int i = 0; i <= 7; ++i) { A4(i, i) = 9.0; }
	for (int i = 0; i <= 5; ++i) A4(i, i +2) = -2.0;
	for (int i = 2; i <= 7; ++i) A4(i, i - 2) = 4.0;
	for (int i = 0; i <= 4; ++i) A4(i, i + 3) = -2.0;
	for (int i = 3; i <= 7; ++i) A4(i, i - 3) = -1.0;*/
	/*std::tuple<Matrix, Matrix, Matrix> result = lu_row_pivoting(A1);
	Matrix P1(std::get<0>(result)), L1(std::get<1>(result)), U1(std::get<2>(result));
	//WriteToCsv("P1.csv", P1); WriteToCsv("L1.csv", L1); WriteToCsv("U1.csv", U1);
	Vector b1(row);
	for (int i = 0; i <= 8; ++i) b1(i) = std::sqrt(0.5*double(i*i) - double(i) + 6.0);
	Vector v1(lu_row_pivoting_solver(A1, b1));
	Vector diff(b1 - A1*v1);
	//WriteToCsv("V1andResid.csv", v1, Norm2(diff));
	Matrix invL1(inverse_triangle(L1, false)), invU1(inverse_triangle(U1, true)), invA1(inverse_regular(A1));
	Vector v2(invA1*b1);
	Vector diff2(b1 - A1*v2);
	//WriteToCsv("V2andResid.csv", v2, Norm2(diff2));
	//Q2
	Matrix A2(A1.transpose()*A1);
	Matrix U2(Cholesky_matrix(A2));
	//WriteToCsv("U2.csv", U2);
	Vector b2(row);
	for (int i = 0; i <= 8; ++i) b2(i) = double(i*i - 9.0) / double(i + 5.0);
	Vector x2(Linear_chol_solver(A2, b2));
	Vector diff3(b2 - A2*x2);
	//WriteToCsv("x2.csv", x2, Norm2(diff3));
	Matrix A3(A1.transpose() + A1);
	Vector x3(lu_row_pivoting_solver(A3, b2));
	Vector diff4(b2 - A3*x3);
	//WriteToCsv("x3.csv", x3, Norm2(diff4));
	WriteToCsv("invA1.csv", invA1); WriteToCsv("invL1.csv", invL1); WriteToCsv("invU1.csv", invU1);
	*/
	/*Vector b4(row);
	for (int i = 0; i <= 7; ++i) b4(i) = double(4.0*i - 3.0) / double(2 * i*i + 1.0);
	Vector x0(row); x0.setZero();
	std::tuple<Vector, int> jaco(Jacobi_residual(A4, b4, tol, x0)), gauss(GaussSiedel_residual(A4, b4, tol, x0)), sor1(SOR_residual(A4, b4, tol, 0.95, x0)), sor2(SOR_residual(A4, b4, tol, 1.21, x0));
	WriteToCsv("jaco_result.csv", std::get<0>(jaco), std::get<1>(jaco));
	WriteToCsv("gauss_result.csv", std::get<0>(gauss), std::get<1>(gauss));
	WriteToCsv("sor1_result.csv", std::get<0>(sor1), std::get<1>(sor1));
	WriteToCsv("sor2_result.csv", std::get<0>(sor2), std::get<1>(sor2));*/
	/*std::ifstream file("data.csv");
	IniMatrixCsv(file, Price);
	for (int i = 0; i < row; ++i) {
		for (int j = 0; j < col; ++j) {
			X(i, j) = std::log(Price(i + 1, j) / Price(i, j));
		}
	}
	//std::cout << X << "\n";
	
	Matrix Cov(GetCov(X));
	Matrix Ux(Cholesky_matrix(Cov));
	Vector y(X.col(5));
	Matrix X_ind(row,col - 1);
	int count = 0;
	for (int j = 0; j < col; ++j) {
		if (j != 5) { X_ind.col(count++) = X.col(j); }
	}
	Vector Coef(LinearRegressionCoef(y, X_ind));
	Vector diff(y - X_ind*Coef);
	WriteToCsv("noConstCoef.csv", Coef,Norm2(diff));
	Matrix X_ind2(row,col);
	X_ind2.col(0).setOnes();
	for (int j = 0; j < col - 1; ++j) X_ind2.col(j + 1) = X_ind.col(j);
	Vector Coef2(LinearRegressionCoef(y, X_ind2));
	Vector diff2(y - X_ind2*Coef2);
	WriteToCsv("ConstCoef.csv", Coef2, Norm2(diff2));
	Matrix M(Ux.transpose());
	WriteToCsv("M.csv", M);*/
	return 0;
}
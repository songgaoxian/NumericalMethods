//implement efficient cubic spline
#ifndef Cubic_Spline_HPP
#define Cubic_Spline_HPP
#include"LinearSolver.hpp"
//matrix中的所有时间为多少个月
Matrix DiscountToRate(Matrix discount, double overnight) {
	int row = discount.rows() + 1, col = 2;
	Matrix Rate(row, col);
	Rate(0, 0) = 0.0;
	Rate(0, 1) = overnight;
	for (int i = 1; i < row; ++i) {
		Rate(i, 0) = double(discount(i - 1, 0))/12.0;//time
		Rate(i,1) = -std::log(discount(i - 1, 1)) / Rate(i, 0);
	}
	std::cout << "Rate:\n" << Rate << "\n";
	return Rate;
}
std::tuple<BandMatrix,Vector> GetTridiagonalSys(Matrix Rate) {
	Vector X(Rate.col(0)), V(Rate.col(1));
	int row = Rate.rows() - 2;
	BandMatrix M(row, row, 1, 1);
	Vector Right(row);
	for (int i = 1; i <= row; ++i) {
		item(M,i - 1, i - 1) = (X(i + 1) - X(i - 1))*2.0 ;
		Right(i - 1) = ((V(i + 1) - V(i)) / (X(i + 1) - X(i)) - (V(i) - V(i - 1)) / (X(i) - X(i - 1)))*6.0;
	}
	for (int i = 2; i <= row; ++i) {
		item(M, i - 1, i - 2) = (X(i) - X(i - 1));
	}
	for (int i = 1; i <= row - 1; ++i) {
		item(M, i - 1, i) = (X(i + 1) - X(i));
	}
	std::cout << "M="<<M.toDenseMatrix()<< "\n";
	std::tuple<BandMatrix, Vector> Result(M, Right);
	return Result;
}
Matrix Cubic_Time_Rate_Discount(Matrix Rate) {
	std::tuple<BandMatrix, Vector> trisys(GetTridiagonalSys(Rate));
	BandMatrix M(std::get<0>(trisys));
	Vector Right(std::get<1>(trisys));
	Vector W_temp(Linear_tri_solver(M, Right));
	int row = Rate.rows();
	Vector W(row), X(Rate.col(0)), V(Rate.col(1)), A(row - 1), B(row - 1), C(row - 1), D(row - 1);
	W(0) = 0; W(row - 1) = 0;
	for (int i = 1; i < row - 1; ++i) W(i) = W_temp(i - 1);
	for (int i = 1; i <= row - 1; ++i) {
		C(i - 1) = (W(i - 1)*X(i) - W(i)*X(i - 1)) / 2.0 / (X(i) - X(i - 1));
		D(i - 1) = (W(i) - W(i - 1)) / 6.0 / (X(i) - X(i - 1));
		double q = V(i - 1) - C(i - 1)*X(i - 1)*X(i - 1) - D(i - 1)*std::pow(X(i - 1), 3);
		double r = V(i) - C(i - 1)*X(i)*X(i) - D(i - 1)*std::pow(X(i), 3);
		A(i - 1) = (q*X(i) - r*X(i - 1)) / (X(i) - X(i - 1));
		B(i - 1) = (r - q) / (X(i) - X(i - 1));
	}
	int max = std::round(12.0*X(row - 1));//max months
	Matrix Result(max, 4);//store zero curve
	Result(max - 1, 0) = X(row - 1);
	Result(max - 1, 1) = Rate(row - 1);
	Result(max - 1, 2) = std::exp(-Rate(row - 1)*X(row - 1));
	Result(max - 1, 3) = max;
	for (int i = 1; i < max; ++i) {
		double c_time = double(i) / 12.0;
		for (int j = 1; j <= row - 1; ++j) {
			if (c_time <= X(j)) {
				double rate = A(j - 1) + B(j - 1)*c_time + C(j - 1)*c_time*c_time + D(j - 1)*c_time*c_time*c_time;
				Result(i - 1, 0) = c_time;
				Result(i - 1, 1) = rate;
				Result(i - 1, 2) = std::exp(-rate*c_time);
				Result(i - 1, 3) = i;
				break;
			}
		}
	}
	return Result;
}
double Price(Matrix Time_CashFlow, Matrix Rate) {
	Matrix ZeroCurve(Cubic_Time_Rate_Discount(Rate));
	int total_month = ZeroCurve.rows();
	int num_cf = Time_CashFlow.rows();
	double price = 0;
	for (int i = 0; i < num_cf; ++i) {
		int current_month = std::round(Time_CashFlow(i, 0));
		double cf = Time_CashFlow(i, 1);
		for (int j = 1; j <= total_month; ++j) {
			if (current_month == j) {
				price += cf*ZeroCurve(j - 1, 2);
				break;
			}
		}
	}
	return price;
}



#endif
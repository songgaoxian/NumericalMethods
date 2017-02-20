//implement classes to solve HeatEquation
#ifndef FiniteDifference_HPP
#define FiniteDifference_HPP
#include<functional>
#include<cmath>
#include "C:\Users\SAMSUNG\Desktop\MTH9821\LU_Forward_Backward\LU_Forward_Backward\Iteration.hpp"
#include "C:\Users\SAMSUNG\Desktop\MTH9821\LU_Forward_Backward\LU_Forward_Backward\LU_Solver.hpp"

std::function<double(double x)> fx1 = [](double x) {return std::exp(x); };
std::function<double(double t)> left1 = [](double t) {return std::exp(t - 2.0); };
std::function<double(double t)> right1 = [](double t) {return std::exp(t + 2.0); };
std::function<double(double, double)> uexact = [](double t, double x) {return std::exp(t + x); };

using Funct = std::function<double(double)>;

double cdfNorm(double x) {
	return 0.5 + 0.5*std::erf(x / std::sqrt(2.0));
}

std::tuple<double, double, double, double> BSPricing(double s, double sigma, double q, double r, double k, double t, bool isCall) {
	double d1 = (std::log(s / k) + (r - q + 0.5*sigma*sigma)*t) / sigma / std::sqrt(t);
	double d2 = d1 - sigma*std::sqrt(t);
	double price = isCall ? s*std::exp(-q*t)*cdfNorm(d1) - k*std::exp(-r*t)*cdfNorm(d2) : k*std::exp(-r*t)*cdfNorm(-d2) - s*std::exp(-q*t)*cdfNorm(-d1);
	double delta = isCall ? std::exp(-q*t)*cdfNorm(d1) : -std::exp(-q*t)*cdfNorm(-d1);
	double gamma = std::exp(-q*t) / s / sigma / std::sqrt(t) / std::sqrt(2.0 * M_PI)*std::exp(-d1*d1 / 2.0);
	double common = -s*sigma*std::exp(-q*t)*std::exp(-d1*d1 / 2.0) / 2.0 / std::sqrt(2.0*M_PI*t);
	double theta = isCall ? common + q*s*std::exp(-q*t)*cdfNorm(d1) - r*k*std::exp(-r*t)*cdfNorm(d2) : common - q*s*std::exp(-q*t)*cdfNorm(-d1) + r*k*std::exp(-r*t)*cdfNorm(-d2);
	std::tuple<double, double, double, double> result(price, delta, gamma, theta);
	return result;
}
double BSPrice(double s, double sigma, double q, double r, double k, double t, bool isCall) {
	double d1 = (std::log(s / k) + (r - q + 0.5*sigma*sigma)*t) / sigma / std::sqrt(t);
	double d2 = d1 - sigma*std::sqrt(t);
	return isCall ? s*std::exp(-q*t)*cdfNorm(d1) - k*std::exp(-r*t)*cdfNorm(d2) : k*std::exp(-r*t)*cdfNorm(-d2) - s*std::exp(-q*t)*cdfNorm(-d1);
}
Vector VarRedResult(Vector us_approx, Vector euro_approx, Vector euro_exact) {
	Vector VarRed = us_approx - (euro_approx - euro_exact);
	return VarRed;
}

class HeatEquation {
public:
	double xleft, xright, tfinal;
	Funct fx, gleft, gright;
	HeatEquation(double xleft_, double xright_, double tfinal_, const Funct& fx_, const Funct& gleft_, const Funct& gright_) :
		xleft(xleft_), xright(xright_), tfinal(tfinal_), fx(fx_), gleft(gleft_), gright(gright_) {}
	HeatEquation() {}
};

class FiniteDifference{
public:
	HeatEquation heq;
	Matrix Result;
	double a;//it is actually alpha
	double tol, w;
	int M; //nums on t
	int N; //nums on x
	double dt, dx;
	FiniteDifference(int M_, int N_, double xleft_, double xright_, double tfinal_, const Funct& fx_, const Funct& gleft_, const Funct& gright_,double tol_,double w_) :
		M(M_), N(N_), heq(xleft_, xright_, tfinal_, fx_, gleft_, gright_),Result(M_+1,N_+1),tol(tol_), w(w_) {
		dt = heq.tfinal / double(M);
		dx = (heq.xright - heq.xleft) / double(N);
		a = dt / (dx*dx);
		double x, t;
		for (int i = 0; i <= N; ++i) {
			x = dx*i + heq.xleft;
			Result(0, i) = heq.fx(x);
		}
		for (int i = 0; i <= M; ++i) {
			t = dt*i;
			Result(i, 0) = heq.gleft(t);
			Result(i, N) = heq.gright(t);
		}
	}
	FiniteDifference() {}
	int GetM() const { return M; }
	int GetN() const { return N; }
	double GetA() const { return a; }
	Matrix GetResult() const { return Result; }
	
	Matrix ForwardEuler() {
		for(int i=1;i<=M;++i)
		  for (int j = 1; j <= N - 1; ++j) {
			  Result(i, j) = a*Result(i - 1, j - 1) + (1.0 - 2 * a)*Result(i - 1, j) + a*Result(i - 1, j + 1);
			}
		return Result;
	}
	Matrix BackwardEuler(bool isLU) {
		Vector Um(N - 1), Um1(N - 1), zero(N - 1), b(N - 1);
		zero.setZero();
		BandMatrix A(N - 1, N - 1, 1, 1);
		for (int i = 0; i < N - 1; ++i) item(A, i, i) = 1 + 2 * a;
		for (int i = 0; i < N - 2; ++i) item(A, i, i + 1) = -a;
		for (int i = 1; i < N - 1; ++i) item(A, i, i - 1) = -a;
		
		for (int i = 0; i < M; ++i) {
			for (int j = 1; j <= N - 1; ++j) {
				Um(j-1) = Result(i, j);
				if (j == 1) zero(j-1)= a*Result(i + 1, 0);
				if (j == N - 1) zero(j-1) = a*Result(i + 1, N);
			}
			b = Um + zero;
			if (isLU) {
				Um1 = lu_no_pivoting_solver(A, b);
			}
			else {
				auto temp = SOR_consecutive(A, b, tol, w, Um);
				Um1 = std::get<0>(temp);
			}
			for (int j = 1; j <= N - 1; ++j) {
				Result(i + 1, j) = Um1(j - 1);
			}
		}
		return Result;
	}
	Matrix CrankNicolson(bool isLU) {
		BandMatrix A(N - 1, N - 1, 1, 1), B(N - 1, N - 1, 1, 1);
		Vector Um(N - 1), Um1(N - 1), zero(N - 1), b(N - 1);
		zero.setZero();
		for (int i = 0; i < N - 1; ++i) {
			item(A, i, i) = 1 + a;
			item(B, i, i) = 1 - a;
		}
		for (int i = 0; i < N - 2; ++i) {
			item(A, i, i + 1) = -0.5*a;
			item(B, i, i + 1) = 0.5*a;
		}
		for (int i = 1; i < N - 1; ++i) {
			item(A, i, i-1) = -0.5*a;
			item(B, i, i - 1) = 0.5*a;
		}
		for (int i = 0; i < M; ++i) {
			for (int j = 1; j <= N - 1; ++j) {
				Um(j - 1) = Result(i, j);
				if (j == 1) zero(j - 1) = 0.5*a*Result(i + 1, 0)+0.5*a*Result(i,0);
				if (j == N - 1) zero(j - 1) = 0.5*a*Result(i + 1, N)+0.5*a*Result(i,N);
			}
			b = vectorMultBand(B,Um) + zero;
			if (isLU) {
				Um1 = lu_no_pivoting_solver(A, b);
			}
			else {
				auto temp = SOR_consecutive(A, b, tol, w, Um);
				Um1 = std::get<0>(temp);
			}
			for (int j = 1; j <= N - 1; ++j) {
				Result(i + 1, j) = Um1(j - 1);
			}
		}
		return Result;
	}
	Matrix GetExact(std::function<double(double, double)>& exact) {
		Matrix ExactResult(M + 1, N + 1);
		double x, t;
		for(int i=0;i<=M;++i)
			for (int j = 0; j <= N; ++j) {
				t = 0 + dt*i;
				x = heq.xleft + dx*j;
				ExactResult(i, j) = exact(t, x);
			}
		return ExactResult;
	}
	std::tuple<double, double> GetError(std::function<double(double, double)>& exact) {
		Matrix ExactM(GetExact(exact));
		double max_err = 0, error_RMS;
		for (int i = 1; i <= N - 1; ++i)
			max_err = std::max(max_err, std::abs(Result(M, i) - ExactM(M, i)));
		double sum = 0;
		for (int i = 0; i <= N; ++i) sum += std::pow((Result(M, i) - ExactM(M, i)), 2) / std::pow(ExactM(M, i), 2);
		error_RMS = std::sqrt(1.0 / double(N + 1)*sum);
		std::tuple<double, double> error(max_err, error_RMS);
		return error;
	}
};

#endif
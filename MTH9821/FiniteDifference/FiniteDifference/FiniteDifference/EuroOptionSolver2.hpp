#ifndef EuroOptionSolver2_HPP
#define EuroOptionSolver2_HPP

#include "HeatEquationSolver.hpp"
#include<vector>

class EuroOptionSolver2 {
protected:
	double K, sigma, s0, q, T, r, xleft, xright, tfinal, a, b;
	std::function<double(double)> payoff, rightf, leftf;
	double numStd;
	FiniteDifference fdm;
	Matrix U;
	Matrix Result;
	bool isCall;
public:
	EuroOptionSolver2(bool isCall_,double numStd_, double K_, double sigma_, double s0_, double q_, double T_, double r_, int M_, double alpha_, double tol_, double w_) :isCall(isCall_), K(K_), sigma(sigma_), s0(s0_), q(q_), T(T_), r(r_) {
		numStd = numStd_;
		a = (r - q) / (sigma*sigma) - 0.5;
		b = std::pow(((r - q) / (sigma*sigma) + 0.5), 2) + 2 * q / (sigma*sigma);
		xleft = std::log(s0 / K) + (r - q - 0.5*sigma*sigma)*T - numStd * sigma*std::sqrt(T);
		xright = std::log(s0 / K) + (r - q - 0.5*sigma*sigma)*T + numStd * sigma*std::sqrt(T);
		tfinal = 0.5*T*sigma*sigma;
		double dt = tfinal / double(M_);
		int N_ = int((xright - xleft) / std::sqrt(dt / alpha_));
		//need to be customized for different option type
		payoff = [=](double x0) {
			if (isCall == false)
				return K*std::exp(a*x0)*std::max(1 - std::exp(x0), 0.0);
			else
				return K*std::exp(a*x0)*std::max(std::exp(x0) - 1, 0.0);
		};
		rightf = [=](double t0) {
			if (isCall == false)
				return 0.0;
			else
				return K*std::exp(a*xright + b*t0)*(-std::exp(-2 * r*t0 / (sigma*sigma)) + std::exp(xright - 2 * q*t0 / sigma / sigma));
		};
		leftf = [=](double t0) {
			if (isCall == false)
				return K*std::exp(a*xleft + b*t0)*(std::exp(-2 * r*t0 / (sigma*sigma)) - std::exp(xleft - 2 * q*t0 / sigma / sigma));
			else
				return 0.0;
		};
		fdm = FiniteDifference(M_, N_, xleft, xright, tfinal, payoff, leftf, rightf, tol_, w_);
	}
	std::tuple<double, double, double, double, double, double> Domain() {
		std::tuple<double, double, double, double, double, double> result(fdm.a, xleft, xright, fdm.N, fdm.dx, fdm.dt);
		return result;
	}
	//need to be adjusted if the option is call
	double BS_ExactPrice(double s1) {
		double d1 = (std::log(s1 / K) + (r - q + 0.5*sigma*sigma)*T) / (sigma*std::sqrt(T));
		double d2 = d1 - sigma*std::sqrt(T);
		if (isCall == false)
			return K*std::exp(-r*T)*cdfNorm(-d2) - s1*std::exp(-q*T)*cdfNorm(-d1);
		else
			return s1*std::exp(-q*T)*cdfNorm(d1) - K*std::exp(-r*T)*cdfNorm(d2);
	}
	void ProcessResult() {
		double xj = xleft, tau = 0;
		Result = U;
		for (int i = 0; i <= fdm.M; ++i)
			for (int j = 0; j <= fdm.N; ++j) {
				xj = xleft + fdm.dx*j;
				tau = fdm.dt*i;
				Result(i, j) = std::exp(-a*xj - b*tau)*Result(i, j);
			}
	}
	Matrix GetU() const { return U; }
	Matrix ForwardEuler() {
		U = fdm.ForwardEuler();
		ProcessResult();
		return Result;
	}
	Matrix BackwardEuler(bool isLU) {
		U = fdm.BackwardEuler(isLU);
		ProcessResult();
		return Result;
	}
	Matrix CrankNicolson(bool isLU) {
		U = fdm.CrankNicolson(isLU);
		ProcessResult();
		return Result;
	}
	int IndexApprox() {
		int indexi = 0;
		double xi = xleft + indexi*fdm.dx;
		while (xi < std::log(s0 / K)) {
			++indexi;
			xi = xleft + indexi*fdm.dx;
		}
		--indexi;
		return indexi;
	}

	std::tuple<double, double> Approximate0() {
		int indexi = IndexApprox();
		double xi = xleft + indexi*fdm.dx;
		double xiplus = xi + fdm.dx;
		double Vi = Result(fdm.M, indexi);
		double Viplus = Result(fdm.M, indexi + 1);
		double si = K*std::exp(xi);
		double siplus = K*std::exp(xiplus);
		double Vapprox0 = ((siplus - s0)*Vi + (s0 - si)*Viplus) / (siplus - si);
		double error1 = std::abs(Vapprox0 - BS_ExactPrice(s0));
		std::tuple<double, double> approx(Vapprox0, error1);
		return approx;
	}
	std::tuple<double, double> GetUs() {
		int indexi = IndexApprox();
		double Ui = U(fdm.M, indexi);
		double Uiplus = U(fdm.M, indexi + 1);
		std::tuple<double, double> approx(Ui, Uiplus);
		return approx;
	}
	std::tuple<double, double> Approximate1() {
		double xcompute = std::log(s0 / K);
		int indexi = IndexApprox();
		double xi = xleft + indexi*fdm.dx;
		double xiplus = xi + fdm.dx;
		double ucompute = ((xiplus - xcompute)*U(fdm.M, indexi) + (xcompute - xi)*U(fdm.M, indexi + 1)) / (xiplus - xi);
		double Vapprox1 = std::exp(-a*xcompute - b*tfinal)*ucompute;
		double error2 = std::abs(Vapprox1 - BS_ExactPrice(s0));
		std::tuple<double, double> approx(Vapprox1, error2);
		return approx;
	}
	double EndError() {
		int count = 0;
		double sum = 0;
		for (int j = 0; j <= fdm.N; ++j) {
			double xj = xleft + j*fdm.dx;
			double Vapprox = Result(fdm.M, j);
			double realval = BS_ExactPrice(K*std::exp(xj));
			if (realval > 0.00001*s0) {
				++count;
				sum += std::pow((Vapprox - realval) / realval, 2);
			}
		}
		double error3 = std::sqrt(sum / double(count));
		return error3;
	}
	std::tuple<double, double, double> DeltaGammaTheta() {
		int indexi = IndexApprox();
		double xi = xleft + indexi*fdm.dx;
		double xiplus = xi + fdm.dx;
		double Vi = Result(fdm.M, indexi);
		double Viplus = Result(fdm.M, indexi + 1);
		double si = K*std::exp(xi);
		double siplus = K*std::exp(xiplus);
		double delta = (Viplus - Vi) / (siplus - si);
		double xiplus2 = xi + 2 * fdm.dx;
		double Viplus2 = Result(fdm.M, indexi + 2);
		double ximinus = xi - fdm.dx;
		double Viminus = Result(fdm.M, indexi - 1);
		double siplus2 = K*std::exp(xiplus2);
		double siminus = K*std::exp(ximinus);
		double diff1 = (Viplus2 - Viplus) / (siplus2 - siplus);
		double diff2 = (Vi - Viminus) / (si - siminus);
		double deno1 = 0.5*(siplus2 + siplus);
		double deno2 = 0.5*(si + siminus);
		double gamma = (diff1 - diff2) / (deno1 - deno2);
		double dtt = 2 * fdm.dt / (sigma*sigma);
		double Vidtt = Result(fdm.M - 1, indexi);
		double Viplusdtt = Result(fdm.M - 1, indexi + 1);
		double Vapproxdtt = ((siplus - s0)*Vidtt + (s0 - si)*Viplusdtt) / (siplus - si);
		double Vapprox0 = ((siplus - s0)*Vi + (s0 - si)*Viplus) / (siplus - si);
		double theta = -(Vapprox0 - Vapproxdtt) / dtt;
		std::tuple<double, double, double> greeks(delta, gamma, theta);
		return greeks;
	}
	std::tuple<double, double, double> ExactGreeks() {
		auto tup = BSPricing(s0, sigma, q, r, K, T, isCall);
		std::tuple<double, double, double> result(std::get<1>(tup), std::get<2>(tup), std::get<3>(tup));
		return result;
	}
	std::vector<Vector> GetVs() {
		std::vector<Vector> result0;
		Vector v1(4), v2(2);
		int indexi = IndexApprox();
		double xi = xleft + indexi*fdm.dx;
		double xiplus = xi + fdm.dx;
		double Vi = Result(fdm.M, indexi);
		double Viplus = Result(fdm.M, indexi + 1);
		double si = K*std::exp(xi);
		double siplus = K*std::exp(xiplus);
		double delta = (Viplus - Vi) / (siplus - si);
		double xiplus2 = xi + 2 * fdm.dx;
		double Viplus2 = Result(fdm.M, indexi + 2);
		double ximinus = xi - fdm.dx;
		double Viminus = Result(fdm.M, indexi - 1);
		v1(0) = Viplus2;
		v1(1) = Viplus;
		v1(2) = Vi;
		v1(3) = Viminus;
		double siplus2 = K*std::exp(xiplus2);
		double siminus = K*std::exp(ximinus);
		double diff1 = (Viplus2 - Viplus) / (siplus2 - siplus);
		double diff2 = (Vi - Viminus) / (si - siminus);
		double deno1 = 0.5*(siplus2 + siplus);
		double deno2 = 0.5*(si + siminus);
		double gamma = (diff1 - diff2) / (deno1 - deno2);
		double dtt = 2 * fdm.dt / (sigma*sigma);
		double Vidtt = Result(fdm.M - 1, indexi);
		double Viplusdtt = Result(fdm.M - 1, indexi + 1);
		v2(0) = Viplusdtt;
		v2(1) = Vidtt;
		result0.push_back(v1);
		result0.push_back(v2);
		return result0;
	}
	std::vector<Vector> GetExactVs() {
		std::vector<Vector> result0;
		Vector v1(4), v2(2);
		int indexi = IndexApprox();
		double xi = xleft + indexi*fdm.dx;
		double xiplus = xi + fdm.dx;
		double Vi = Result(fdm.M, indexi);
		double Viplus = Result(fdm.M, indexi + 1);
		double si = K*std::exp(xi);
		double siplus = K*std::exp(xiplus);
		double delta = (Viplus - Vi) / (siplus - si);
		double xiplus2 = xi + 2 * fdm.dx;
		double Viplus2 = Result(fdm.M, indexi + 2);
		double ximinus = xi - fdm.dx;
		double Viminus = Result(fdm.M, indexi - 1);
		double siplus2 = K*std::exp(xiplus2);
		double siminus = K*std::exp(ximinus);
		v1(0) = BSPrice(siplus2, sigma, q, r, K, T, isCall);
		v1(1) = BSPrice(siplus, sigma, q, r, K, T, isCall);
		v1(2) = BSPrice(si, sigma, q, r, K, T, isCall);
		v1(3) = BSPrice(siminus, sigma, q, r, K, T, isCall);
		double dtt = 2 * fdm.dt / (sigma*sigma);
		double newT = T - dtt;
		v2(0) = BSPrice(siplus, sigma, q, r, K, newT, isCall);
		v2(1) = BSPrice(si, sigma, q, r, K, newT, isCall);
		result0.push_back(v1);
		result0.push_back(v2);
		return result0;
	}
	std::tuple<double, double, double, double, double, double, double, double> Summary() {
		std::tuple<double, double, double, double, double, double, double, double>
			summary(std::get<0>(Approximate0()), std::get<1>(Approximate0()), std::get<0>(Approximate1()), std::get<1>(Approximate1()), EndError(),
				std::get<0>(DeltaGammaTheta()), std::get<1>(DeltaGammaTheta()), std::get<2>(DeltaGammaTheta()));
		return summary;
	}
	std::tuple<double, double, double, double> Summary2() {
		std::tuple<double, double, double, double>
			summary(std::get<0>(Approximate0()),
				std::get<0>(DeltaGammaTheta()), std::get<1>(DeltaGammaTheta()), std::get<2>(DeltaGammaTheta()));
		return summary;
	}
};
#endif

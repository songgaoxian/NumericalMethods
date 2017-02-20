#ifndef AmericanOptionSolver2_HPP
#define AmericanOptionSolver2_HPP

#include "EuroOptionSolver.hpp"
#include"EuroOptionSolver2.hpp"

class AmericanOptionSolver2 {
private:
	double K, sigma, s0, q, T, r, xleft, xright, tfinal, a, b, alpha, xcompute1;
	double numStd;
	std::function<double(double)> payoff, rightf, leftf;
	double tol, w;
	Matrix U;
	Matrix Result;
	bool isCall;
	int M, N;
	double dt, dx;
	double exactP;
	FiniteDifference fdm;
public:
	AmericanOptionSolver2(bool isCall_,double numStd_, double exactP_, double K_, double sigma_, double s0_, double q_, double T_, double r_, int M_, double alpha_, double tol_, double w_) :
		isCall(isCall_),numStd(numStd_), exactP(exactP_), M(M_), K(K_), sigma(sigma_), s0(s0_), q(q_), T(T_), r(r_), alpha(alpha_), tol(tol_), w(w_) {
		a = (r - q) / (sigma*sigma) - 0.5;
		xcompute1 = std::log(s0 / K);
		b = std::pow(((r - q) / (sigma*sigma) + 0.5), 2) + 2 * q / (sigma*sigma);
		xleft = std::log(s0 / K) + (r - q - 0.5*sigma*sigma)*T - numStd * sigma*std::sqrt(T);
		xright = std::log(s0 / K) + (r - q - 0.5*sigma*sigma)*T + numStd * sigma*std::sqrt(T);
		tfinal = 0.5*T*sigma*sigma;
		dt = tfinal / double(M_);
		N = int((xright - xleft) / std::sqrt(dt / alpha));
		dx = (xright - xleft) / double(N);
		alpha = dt / (dx*dx);
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
				return K*std::exp(a*xright + b*t0)*(std::exp(xright) - 1.0);
		};
		leftf = [=](double t0) {
			if (isCall == false)
				return K*std::exp(a*xleft + b*t0)*(1.0 - std::exp(xleft));
			else
				return 0.0;
		};
		fdm = FiniteDifference(M, N, xleft, xright, tfinal, payoff, leftf, rightf, tol_, w_);
		U = fdm.GetResult();
		Result = U;
	}
	std::tuple<double, double, double, double, double, double, double, double> DomainInfo() {
		std::tuple<double, double, double, double, double, double, double, double> result(fdm.a, N, xleft, xright, xcompute1, tfinal, dt, dx);
		return result;
	}
	double Early_Ex_Premium(double x, double tau) {
		if (isCall == false)
			return K*std::exp(a*x + b*tau)*std::max(1.0 - std::exp(x), 0.0);
		else
			return K*std::exp(a*x + b*tau)*std::max(std::exp(x) - 1.0, 0.0);
	}
	void ProcessResult() {
		double xj = xleft, tau = 0;
		for (int i = 0; i <= M; ++i)
			for (int j = 0; j <= N; ++j) {
				xj = xleft + dx*j;
				tau = dt*i;
				Result(i, j) = std::exp(-a*xj - b*tau)*U(i, j);
			}
	}
	Matrix ForwardEuler() {
		double x, tau;
		for (int i = 1; i <= M; ++i)
			for (int j = 1; j <= N - 1; ++j) {
				tau = dt*i;
				x = xleft + dx*j;
				U(i, j) = alpha*U(i - 1, j - 1) + (1.0 - 2 * alpha)*U(i - 1, j) + alpha*U(i - 1, j + 1);
				U(i, j) = std::max(U(i, j), Early_Ex_Premium(x, tau));
			}
		ProcessResult();
		return Result;
	}

	Matrix CrankNicolson() {
		BandMatrix B(N - 1, N - 1, 1, 1), A(N - 1, N - 1, 1, 1);
		Vector Um(N - 1), Um1(N - 1), zero(N - 1), b(N - 1);
		zero.setZero();
		for (int i = 0; i < N - 1; ++i) {
			item(A, i, i) = 1 + alpha;
			item(B, i, i) = 1 - alpha;
		}
		for (int i = 0; i < N - 2; ++i) {
			item(A, i, i + 1) = -0.5*alpha;
			item(B, i, i + 1) = 0.5*alpha;
		}
		for (int i = 1; i < N - 1; ++i) {
			item(A, i, i - 1) = -0.5*alpha;
			item(B, i, i - 1) = 0.5*alpha;
		}
		for (int i = 0; i < M; ++i) {
			for (int j = 1; j <= N - 1; ++j) {
				Um(j - 1) = U(i, j);
				if (j == 1) zero(j - 1) = 0.5*alpha*U(i + 1, 0) + 0.5*alpha*U(i, 0);
				if (j == N - 1) zero(j - 1) = 0.5*alpha*U(i + 1, N) + 0.5*alpha*U(i, N);
			}
			Vector diff;
			b = vectorMultBand(B, Um) + zero;
			do {
				double origin = U(i + 1, 0);
				for (int j = 1; j <= N - 1; ++j) {
					double tau = dt*i + dt;
					double x = xleft + dx*j;
					if (j == 1)
						Um1(j - 1) = (1.0 - w)*Um(j - 1) + w*alpha / (2.0 * (1.0 + alpha))*(Um(j)) + w / (1.0 + alpha)*b(j - 1);
					else if (j == N - 1)
						Um1(j - 1) = (1.0 - w)*Um(j - 1) + w*alpha / (2.0 * (1.0 + alpha))*(Um1(j - 2)) + w / (1.0 + alpha)*b(j - 1);
					else
						Um1(j - 1) = (1.0 - w)*Um(j - 1) + w*alpha / (2.0 * (1.0 + alpha))*(Um1(j - 2) + Um(j)) + w / (1.0 + alpha)*b(j - 1);
					Um1(j - 1) = std::max(Um1(j - 1), Early_Ex_Premium(x, tau));
				}
				diff = Um1 - Um;
				Um = Um1;
			} while (Norm2(diff) >= tol);
			for (int j = 1; j <= N - 1; ++j) {
				U(i + 1, j) = Um1(j - 1);
			}
		}
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
		double error1 = std::abs(Vapprox0 - exactP);
		std::tuple<double, double> approx(Vapprox0, error1);
		return approx;
	}
	std::tuple<double, double> GetUs() {
		int indexi = IndexApprox();
		double xi = xleft + indexi*fdm.dx;
		double xiplus = xi + fdm.dx;
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
		double error2 = std::abs(Vapprox1 - exactP);
		std::tuple<double, double> approx(Vapprox1, error2);
		return approx;
	}
	std::tuple<double, double> ApproximateVarRed(bool useFE) {
		EuroOptionSolver2 euroSolver(isCall,numStd, K, sigma, s0, q, T, r, M, alpha, tol, w);
		double PvarRed, error;
		if (useFE) {
			ForwardEuler();
			euroSolver.ForwardEuler();
			PvarRed = std::get<0>(Approximate0()) + euroSolver.BS_ExactPrice(s0) - std::get<0>(euroSolver.Approximate0());
			error = std::abs(PvarRed - exactP);
		}
		else {
			CrankNicolson();
			euroSolver.CrankNicolson(true);
			PvarRed = std::get<0>(Approximate0()) + euroSolver.BS_ExactPrice(s0) - std::get<0>(euroSolver.Approximate0());
			error = std::abs(PvarRed - exactP);
		}
		std::tuple<double, double> approx(PvarRed, error);
		return approx;
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
		double dtt = 2 * fdm.dt / (sigma*sigma);
		double Vidtt = Result(fdm.M - 1, indexi);
		double Viplusdtt = Result(fdm.M - 1, indexi + 1);
		v2(0) = Viplusdtt;
		v2(1) = Vidtt;
		result0.push_back(v1);
		result0.push_back(v2);
		return result0;
	}
	std::tuple<Vector, double> GetSsAndDtt() {
		Vector v1(4);
		int indexi = IndexApprox();
		double xi = xleft + indexi*fdm.dx;
		double xiplus = xi + fdm.dx;
		double Vi = Result(fdm.M, indexi);
		double Viplus = Result(fdm.M, indexi + 1);
		double si = K*std::exp(xi);
		double siplus = K*std::exp(xiplus);
		double xiplus2 = xi + 2 * fdm.dx;
		double Viplus2 = Result(fdm.M, indexi + 2);
		double ximinus = xi - fdm.dx;
		double Viminus = Result(fdm.M, indexi - 1);
		double siplus2 = K*std::exp(xiplus2);
		double siminus = K*std::exp(ximinus);
		v1(0) = siplus2;
		v1(1) = siplus;
		v1(2) = si;
		v1(3) = siminus;
		double dtt = 2 * fdm.dt / (sigma*sigma);
		double Vidtt = Result(fdm.M - 1, indexi);
		double Viplusdtt = Result(fdm.M - 1, indexi + 1);
		std::tuple<Vector, double> result0(v1, dtt);
		return result0;
	}
	std::vector<Vector> GetExactVs() {
		auto tup = GetSsAndDtt();
		double dtt = std::get<1>(tup);
		Vector Ss = std::get<0>(tup);
		Vector v1(4), v2(2);
		v1(0) = BSPrice(Ss(0), sigma, q, r, K, T, isCall);
		v1(1) = BSPrice(Ss(1), sigma, q, r, K, T, isCall);
		v1(2) = BSPrice(Ss(2), sigma, q, r, K, T, isCall);
		v1(3) = BSPrice(Ss(3), sigma, q, r, K, T, isCall);
		v2(0) = BSPrice(Ss(1), sigma, q, r, K, T - dtt, isCall);
		v2(1) = BSPrice(Ss(2), sigma, q, r, K, T - dtt, isCall);
		std::vector<Vector> result0;
		result0.push_back(v1);
		result0.push_back(v2);
		return result0;
	}
	int GetNopt(int m) {
		double tau = m*dt;
		CrankNicolson();
		int Nopt = 0;
		for (int i = 0; i <= N; ++i) {
			double x = xleft + dx*i;
			if (U(m, i) > Early_Ex_Premium(x, tau)) {
				Nopt = i - 1;
				return Nopt;
			}
		}
		if (m == 0) {
			for (int i = 0; i <= N; ++i) {
				double x = xleft + dx*i;
				if (U(m, i) == 0) {
					Nopt = i - 1;
					return Nopt;
				}
			}
		}
		return Nopt;
	}
	Vector Early_Exercise_Domain() {
		Vector Domain(M + 1);
		CrankNicolson();
		for (int i = 0; i <= M; ++i) {
			int index = GetNopt(i);
			double x = xleft + dx*index;
			double xplus = x + dx;
			Domain(i) = 0.5*(K*std::exp(x) + K*std::exp(xplus));
		}
		return Domain;
	}
	Matrix GetU() const { return U; }
	std::tuple<double, double, double, double, double, double, double, double, double> Summary(bool useFE) {
		std::tuple<double, double, double, double, double, double, double, double, double>
			summary(std::get<0>(Approximate0()), std::get<1>(Approximate0()), std::get<0>(Approximate1()), std::get<1>(Approximate1()), std::get<0>(ApproximateVarRed(useFE)),
				std::get<1>(ApproximateVarRed(useFE)), std::get<0>(DeltaGammaTheta()), std::get<1>(DeltaGammaTheta()), std::get<2>(DeltaGammaTheta()));
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

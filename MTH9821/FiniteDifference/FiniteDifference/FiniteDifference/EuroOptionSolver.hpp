#ifndef EuroOptionSolver_HPP
#define EuroOptionSolver_HPP

#include "HeatEquationSolver.hpp"
#include<vector>

class EuroOptionDividend;
class EuroOptionSolver {
	friend class EuroOptionDividend;
protected:
	double K, sigma, s0, q, T, r, xleft, xright, tfinal, a, b;
	std::function<double(double)> payoff, rightf, leftf;
	FiniteDifference fdm;
	Matrix U;
	Matrix Result;
	bool isCall;
public:
	EuroOptionSolver(bool isCall_,double K_, double sigma_, double s0_, double q_, double T_, double r_, int M_,double alpha_,double tol_, double w_) :isCall(isCall_), K(K_), sigma(sigma_), s0(s0_), q(q_), T(T_), r(r_) {
		a = (r - q) / (sigma*sigma) - 0.5;
		b = std::pow(((r - q) / (sigma*sigma) + 0.5), 2) + 2 * q / (sigma*sigma);
		xleft = std::log(s0 / K) + (r - q - 0.5*sigma*sigma)*T - 3 * sigma*std::sqrt(T);
		xright = std::log(s0 / K) + (r - q - 0.5*sigma*sigma)*T + 3 * sigma*std::sqrt(T);
		tfinal = 0.5*T*sigma*sigma;
		double dt = tfinal / double(M_);
		int N_ = int((xright - xleft) / std::sqrt(dt / alpha_));
		//need to be customized for different option type
		payoff= [=](double x0) {
			if (isCall == false)
				return K*std::exp(a*x0)*std::max(1 - std::exp(x0), 0.0);
			else
				return K*std::exp(a*x0)*std::max(std::exp(x0) - 1, 0.0);
		};
		rightf= [=](double t0) {
			if(isCall==false)
			  return 0.0;
			else
			  return K*std::exp(a*xright + b*t0)*(-std::exp(-2 * r*t0 / (sigma*sigma)) + std::exp(xright - 2 * q*t0 / sigma / sigma));
		};
		leftf=[=](double t0) {
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
		for(int i=0;i<=fdm.M;++i)
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
		v2(0) = BSPrice(siplus,sigma,q,r,K,newT,isCall);
		v2(1) =BSPrice(si,sigma,q,r,K,newT,isCall);
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
class EuroOptionDividend {
protected:
	double t_div, q_div, tau_div, xcompute1, xcompute2, alpha1, dtau1, dx, dtau2;
	int Nleft, Nright, N, M1, M2;
	EuroOptionSolver euro_solver, euro_solver2;
public:
	EuroOptionDividend(double t_div_, double q_div_, bool isCall_, double K_, double sigma_, double s0_, double q_, double T_, double r_, int M_, double alpha, double tol_, double w_) :
		t_div(t_div_), q_div(q_div_), euro_solver(isCall_, K_, sigma_, s0_, q_, T_, r_, M_, alpha, tol_, w_), euro_solver2(isCall_, K_, sigma_, s0_, q_, T_, r_, M_, alpha, tol_, w_) {
		tau_div = 0.5*euro_solver.sigma*euro_solver.sigma*(euro_solver.T - t_div);
		xcompute1 = std::log(euro_solver.s0 / euro_solver.K) + std::log(1 - q_div);
		xcompute2 = std::log(euro_solver2.s0 / euro_solver2.K);
		alpha1 = alpha;
		dtau1 = tau_div / double(M_);
		M1 = M_;
		dx = std::sqrt(dtau1 / alpha);
		double xleft_temp = std::log(euro_solver.s0 / euro_solver.K) + (euro_solver.r - 0.5*euro_solver.sigma*euro_solver.sigma)*euro_solver.T - 3 * euro_solver.sigma*std::sqrt(euro_solver.T);
		double xright_temp = xleft_temp + 6 * euro_solver.sigma*std::sqrt(euro_solver.T);
		Nleft = std::ceil((xcompute1 - xleft_temp) / dx);
		Nright = std::ceil((xright_temp - xcompute1) / dx);
		N = Nleft + Nright;
		euro_solver.xleft = xcompute1 - Nleft*dx;
		euro_solver.xright = xcompute1 + Nright*dx;
		euro_solver.tfinal = tau_div;
		//need to be customized for different option type
		euro_solver.rightf = [=](double t0) {
			if (euro_solver.isCall == false)
				return 0.0;
			else
				return euro_solver.K*std::exp(euro_solver.a*euro_solver.xright + euro_solver.b*t0)*(-std::exp(-2 * euro_solver.r*t0 / (euro_solver.sigma*euro_solver.sigma)) + std::exp(euro_solver.xright - 2 * euro_solver.q*t0 / euro_solver.sigma / euro_solver.sigma));
		};
		euro_solver.leftf = [=](double t0) {
			if (euro_solver.isCall == false)
				return euro_solver.K*std::exp(euro_solver.a*euro_solver.xleft + euro_solver.b*t0)*(std::exp(-2 * euro_solver.r*t0 / (euro_solver.sigma*euro_solver.sigma)) - std::exp(euro_solver.xleft - 2 * euro_solver.q*t0 / euro_solver.sigma / euro_solver.sigma));
			else
				return 0.0;
		};
		euro_solver.fdm = FiniteDifference(M_, N, euro_solver.xleft, euro_solver.xright, tau_div, euro_solver.payoff, euro_solver.leftf, euro_solver.rightf, tol_, w_);
	}
	Vector DomainDiscretization() {
		Vector Result(11);
		euro_solver2.xleft = euro_solver.xleft - std::log(1 - q_div);
		euro_solver2.xright = euro_solver.xright - std::log(1 - q_div);
		double dtau2temp = alpha1*dx*dx;
		M2 = std::ceil((euro_solver2.tfinal - tau_div) / dtau2temp);
		dtau2 = (euro_solver2.tfinal - tau_div) / double(M2);
		double alpha2 = dtau2 / (dx*dx);
		//need to be customized for different option type
		euro_solver2.rightf = euro_solver.rightf;
		euro_solver2.leftf = euro_solver.leftf;
		double tol_ = euro_solver2.fdm.tol;
		double w_ = euro_solver2.fdm.w;
		euro_solver2.fdm = FiniteDifference(M2, N, euro_solver2.xleft, euro_solver2.xright, euro_solver2.tfinal - tau_div, euro_solver2.payoff, euro_solver2.leftf, euro_solver2.rightf, tol_, w_);
		Result(0) = M2;
		Result(1) = alpha2;
		Result(2) = N;
		Result(3) = euro_solver.xleft;
		Result(4) = euro_solver.xright;
		Result(5) = euro_solver2.xleft;
		Result(6) = euro_solver2.xright;
		Result(7) = tau_div;
		Result(8) = dtau1;
		Result(9) = dtau2;
		Result(10) = dx;
		return Result;
	}
	Matrix GetU() const { return euro_solver2.GetU(); }
	Matrix GetU0() const { return euro_solver.GetU(); }
	Matrix ForwardEuler() {
		euro_solver.ForwardEuler();
		auto U1 = euro_solver.GetU();
		DomainDiscretization();
		for (int i = 0; i <= N; ++i) {
			euro_solver2.fdm.Result(0, i) = U1(M1,i);
		}
		return euro_solver2.ForwardEuler();
	}
	Matrix BackwardEuler(bool isLU) {
		euro_solver.BackwardEuler(isLU);
		auto U1 = euro_solver.GetU();
		DomainDiscretization();
		for (int i = 0; i <= N; ++i) {
			euro_solver2.fdm.Result(0, i) = U1(M1, i);
		}
		return euro_solver2.BackwardEuler(isLU);
	}
	Matrix CrankNicolson(bool isLU) {
		euro_solver.CrankNicolson(isLU);
		auto U1 = euro_solver.GetU();
		DomainDiscretization();
		for (int i = 0; i <= N; ++i) {
			euro_solver2.fdm.Result(0, i) = U1(M1, i);
		}
		return euro_solver2.CrankNicolson(isLU);
	}
	Vector ValueGreeks() {
		double Vs0 = euro_solver2.Result(M2, Nleft);
		double Sminus = euro_solver2.K*std::exp(xcompute2 - dx);
		double S0 = euro_solver2.s0;
		double Splus = euro_solver2.K*std::exp(xcompute2 + dx);
		double Vminus = euro_solver2.Result(M2, Nleft - 1);
		double Vplus = euro_solver2.Result(M2, Nleft + 1);
		double Delta = (Vplus - Vminus) / (Splus - Sminus);
		double Gamma = (S0 - Sminus)*Vplus - (Splus - Sminus)*Vs0 + (Splus - S0)*Vminus;
		Gamma = Gamma / (S0 - Sminus) / (Splus - S0) / (0.5*(Splus - Sminus));
		double Vs1 = euro_solver2.Result(M2 - 1, Nleft);
		double dtt = 2 * dtau2 / euro_solver2.sigma / euro_solver2.sigma;
		double Theta = (Vs0 - Vs1) / dtt;
		Vector greeks(5);
		greeks(0) = euro_solver2.U(M2, Nleft);
		greeks(1) = Vs0;
		greeks(2) = Delta;
		greeks(3) = Gamma;
		greeks(4)=Theta;
		return greeks;
	}
};


#endif
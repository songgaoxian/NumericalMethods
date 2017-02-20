//implement EuroBarrier Option
#ifndef EuroBarrier_HPP
#define EuroBarrier_HPP

#include "EuroOptionSolver.hpp"
//only for call
class EuroDownOut : public EuroOptionSolver
{
private:
	double B;
	double xcompute;
	int Nleft;
	int Nright;
public:
	EuroDownOut(double B_, bool isCall_, double K_, double sigma_, double s0_, double q_, double T_, double r_, int M_, double alpha, double tol_, double w_) :
		B(B_), EuroOptionSolver(isCall_, K_, sigma_, s0_, q_, T_, r_, M_, alpha, tol_, w_)
	{
		xleft = std::log(B / K);
		xcompute = std::log(s0 / K);
		double dt = tfinal / double(M_);
		double dxtemp = std::sqrt(dt / alpha);
		Nleft = int((xcompute - xleft) / dxtemp);
		double dx = (xcompute - xleft) / double(Nleft);
		double xrighttemp = std::log(s0 / K) + (r - q - 0.5*sigma*sigma)*T + 3 * sigma*std::sqrt(T);
		Nright = std::ceil((xrighttemp - xcompute) / dx);
		int N_ = Nleft + Nright;
		xright = xcompute + Nright*dx;
		fdm = FiniteDifference(M_, N_, xleft, xright, tfinal, payoff, leftf, rightf, tol_, w_);
	}
	std::tuple<double, double, double, double, double, double,double,double> DomainInfo() {
		std::tuple<double, double, double, double, double, double,double,double> result(fdm.a,fdm.N, xleft, xright,xcompute,tfinal, fdm.dx, fdm.dt);
		return result;
	}
	double VExact(double S) {
		return BS_ExactPrice(S) - std::pow(B / S, 2 * a)*BS_ExactPrice(B*B / S);
	}
	std::tuple<double,double,double,double,double,double> GreeksDataApproximate() {
		double Vapprox = Result(fdm.M, Nleft);
		double uapprox = U(fdm.M, Nleft);
		double exact = VExact(s0);
		double error = std::abs(Vapprox - exact);
		double Sminus = K*std::exp(xcompute - fdm.dx);
		double S0 = s0;
		double Splus = K*std::exp(xcompute + fdm.dx);
		double Vminus = Result(fdm.M, Nleft - 1);
		double Vplus = Result(fdm.M, Nleft + 1);
		double delta = (Vplus - Vminus) / (Splus - Sminus);
		double gamma = ((S0 - Sminus)*Vplus - (Splus - Sminus)*Vapprox + (Splus - S0)*Vminus) / (S0 - Sminus) / (Splus - S0) / (0.5*(Splus - Sminus));
		double Vapprox1 = Result(fdm.M - 1, Nleft);
		double dtt = 2 * fdm.dt / (sigma*sigma);
		double theta = (Vapprox - Vapprox1) / dtt;
		std::tuple<double,double, double, double, double, double> result(uapprox,Vapprox, error, delta, gamma, theta);
		return result;
	}
};


#endif
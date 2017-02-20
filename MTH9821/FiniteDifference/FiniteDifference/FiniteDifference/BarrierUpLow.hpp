#ifndef BarrierUpLow_HPP
#define BarrierUpLow_HPP

#include"EuroBarrier.hpp"

class BarrierUpLow : public EuroOptionSolver
{
private:
	double Bup, Blow;
	double xcompute;
	double payUp, payLow;
public:
	BarrierUpLow(double payUp_,double payLow_,double Bup_,double Blow_, bool isCall_, double K_, double sigma_, double s0_, double q_, double T_, double r_, int M_, double alpha, double tol_, double w_) :
		Bup(Bup_),Blow(Blow_),EuroOptionSolver(isCall_, K_, sigma_, s0_, q_, T_, r_, M_, alpha, tol_, w_)
	{
		payUp = payUp_;
		payLow = payLow_;
		xleft = std::log(Blow / K);
		xcompute = std::log(s0 / K);
		xright = std::log(Bup / K);
		double dt = tfinal / double(M_);
		int N_ = std::floor((xright - xleft) / std::sqrt(dt / alpha));
		rightf = [=](double t0) {
			return payUp*std::exp(a*xright + b*t0);
		};
		leftf = [=](double t0) {
			return payLow*std::exp(a*xleft + b*t0);
		};
		fdm = FiniteDifference(M_, N_, xleft, xright, tfinal, payoff, leftf, rightf, tol_, w_);
	}
	std::tuple<double, double, double, double, double, double, double, double> DomainInfo() {
		std::tuple<double, double, double, double, double, double, double, double> result(fdm.a, fdm.N, xleft, xright, xcompute, tfinal, fdm.dt,fdm.dx);
		return result;
	}
};


#endif
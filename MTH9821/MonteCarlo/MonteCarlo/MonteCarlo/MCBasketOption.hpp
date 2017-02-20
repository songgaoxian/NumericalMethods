//implement BasketOption pricing using Monte Carlo Simulation
#ifndef MCBasketOption_HPP
#define MCBasketOption_HPP

#include "MonteCarlo.hpp"
#include<cmath>

class MCBasketOption : public MonteCarlo
{
protected:
	std::tuple<double, double,double> stock2;
	double rho;
public:
	MCBasketOption(long size_, std::tuple<double, double, double, double, double, double>& stck, int type, std::tuple<double, double,double>& stck2, double corr) :
		MonteCarlo(size_, stck, type), stock2(stck2), rho(corr) {}
	double EuroBasketCall(long N) {
		if (2 * N > normal.size()) throw "normal variables are not enough";
		double S1, S2, V, V_bar = 0;
		double s10 = std::get<0>(stock), s20 = std::get<0>(stock2);
		double k = std::get<1>(stock);
		double t = std::get<2>(stock);
		double sigma1 = std::get<3>(stock), sigma2 = std::get<1>(stock2);
		double q1 = std::get<4>(stock), q2 = std::get<2>(stock2);
		double r = std::get<5>(stock);
		for (int i = 0; i < N; ++i) {
			S1 = s10*std::exp((r - 0.5*sigma1*sigma1 - q1)*t + sigma1*std::sqrt(t)*normal[2 * i]);
			S2 = s20*std::exp((r - 0.5*sigma2*sigma2 - q2)*t + sigma2*std::sqrt(t)*(rho*normal[2 * i] + std::sqrt(1 - rho*rho)*normal[2 * i + 1]));
			V = std::exp(-r*t)*std::max(S1 + S2 - k, 0.0);
			V_bar += V / double(N);
			//if (i % 1000000 == 0) std::cout << i << " ";
		}
		return V_bar;
	}
	double DependentBasketCall(long m, long n) {
		if (2 * n*m > normal.size()) throw "normal variables are not enough";
		double s10 = std::get<0>(stock), s20 = std::get<0>(stock2);
		double k = std::get<1>(stock);
		double t = std::get<2>(stock);
		double sigma1 = std::get<3>(stock), sigma2 = std::get<1>(stock2);
		double q1 = std::get<4>(stock), q2 = std::get<2>(stock2);
		double r = std::get<5>(stock);
		double dt = t / double(m);
		long count = 0;
		double maxsum = s10 + s20;
		double S1, S2, Vi, V_bar = 0;
		for (int i = 0; i < n; ++i) {
			S1 = s10; S2 = s20; maxsum = s10 + s20;
			for (int j = 0; j < m; ++j) {
				double z1 = normal[count];
				++count;
				double z2 = normal[count];
				++count;
				S1 = S1*std::exp((r - q1 - 0.5*sigma1*sigma1)*dt + sigma1*std::sqrt(dt)*z1);
				S2 = S2*std::exp((r - q2 - 0.5*sigma2*sigma2)*dt + sigma2*std::sqrt(dt)*(rho*z1 + std::sqrt(1 - rho*rho)*z2));
				maxsum = std::max(maxsum, S1 + S2);
			}
			Vi = std::max(maxsum - k, 0.0)*std::exp(-r*t);
			V_bar += Vi / double(n);
		}
		return V_bar;
	}

};

class MCBarrierOption : public MonteCarlo
{
protected:
	double Bup, Blow, pay;
public:
	MCBarrierOption(long size_,double Bup_,double Blow_,double pay_, std::tuple<double, double, double, double, double, double>& stck, int type) :
		MonteCarlo(size_, stck, type), Bup(Bup_),Blow(Blow_),pay(pay_) {}

	std::tuple<double,double,double> DependentBarriertCall(long Nk) {
		if (Nk > normal.size()) throw "normal variables are not enough";
		double s10 = std::get<0>(stock);
		double k = std::get<1>(stock);
		double t = std::get<2>(stock);
		long m = std::ceil(std::pow(Nk, 1.0 / 3.0)*std::pow(t, 2.0 / 3.0));
		long n = std::floor(Nk / m);
		double sigma1 = std::get<3>(stock);
		double q1 = std::get<4>(stock);
		double r = std::get<5>(stock);
		double dt = t / double(m);
		long count = 0;
		double S1, Vi, V_bar = 0;
		for (int i = 0; i < n; ++i) {
			S1 = s10;
			bool isValid = true;
			double t_pay;
			for (int j = 0; j < m; ++j) {
				double z1 = normal[count];
				++count;
				S1 = S1*std::exp((r - q1 - 0.5*sigma1*sigma1)*dt + sigma1*std::sqrt(dt)*z1);
				if (isValid && (S1 >= Bup || S1 <= Blow)) { 
					isValid = false;
					t_pay = dt*j + dt;
				}
			}
			if (isValid) {
				Vi = std::max(S1 - k, 0.0)*std::exp(-r*t);
			}
			else {
				Vi = pay*std::exp(-r*t_pay);
			}
			V_bar += Vi / double(n);
		}
		std::tuple<double, double, double> result1(m, n, V_bar);
		return result1;
	}

};

#endif
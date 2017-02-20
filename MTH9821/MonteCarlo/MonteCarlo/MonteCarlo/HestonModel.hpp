//implement HestonModel
#ifndef HestonModel_HPP
#define HestonModel_HPP

#include "MonteCarlo.hpp"

class HestonPricer : public MonteCarlo
{
protected:
	double reversion, mean, std_var, rho, v0;
public:
	HestonPricer(long size_, std::tuple<double, double, double, double, double, double>& stck, int type, double reversion_, double mean_, double stdvar_, double rho_, double v0_) :
		MonteCarlo(size_, stck, type), reversion(reversion_), mean(mean_), std_var(stdvar_), rho(rho_), v0(v0_) {}
	double BSPricerVol(double vol) {
		double s0 = std::get<0>(stock);
		double k = std::get<1>(stock);
		double t = std::get<2>(stock);
		double q = std::get<4>(stock);
		double r = std::get<5>(stock);
		double d1 = (std::log(s0 / k) + (r - q + 0.5*vol*vol)*t) / (vol*std::sqrt(t));
		double d2 = d1 - vol*std::sqrt(t);
		double p = k*std::exp(-r*t)*cdfNorm(-d2) - s0*std::exp(-q*t)*cdfNorm(-d1);
		return p;
	}
	double ImpliedVol(double price,double minvol, double maxvol, double tol) {
		if (BSPricerVol(minvol) > price) throw "reduce volatility lower bound";
		if (BSPricerVol(maxvol) < price) throw "increase volatility upper bound";
		double mid;
		while (maxvol - minvol > tol) {
			mid = 0.5*(maxvol + minvol);
			if (BSPricerVol(mid) > price) maxvol = mid;
			else minvol = mid;
		}
		return mid;
	}
	std::tuple<double, double> HestonPutPricer(long m, long n,double minvol,double maxvol) {
		if (m*n * 2 > normal.size()) throw "normal variables are not enough";
		double s0 = std::get<0>(stock);
		double k = std::get<1>(stock);
		double t = std::get<2>(stock);
		double q = std::get<4>(stock);
		double r = std::get<5>(stock);
		double v, vplus, Vi, V_bar = 0, s;
		double dt = t / double(m);
		long count = -1, zero = 0;
		for (int i = 0; i < n; ++i) {
			s = s0; v = v0; vplus = v0; zero = 0;
			for (int j = 0; j < m; ++j) {
				s = s*std::exp((r - q - 0.5*vplus)*dt + std::sqrt(vplus*dt)*normal[++count]);
				v = vplus - reversion*(vplus - mean*mean)*dt + std_var*std::sqrt(vplus*dt)*(rho*normal[count] + std::sqrt(1 - rho*rho)*normal[++count]);
				vplus = std::max(v, 0.0);
				if (v < 0) zero++;
			}
			Vi = std::exp(-r*t)*std::max(k - s, 0.0);
			V_bar += Vi / double(n);
		}
        double implied_vol = ImpliedVol(V_bar, minvol, maxvol, 1E-9);
		std::tuple<double, double> result(V_bar, implied_vol);
		return result;
	}

};

#endif
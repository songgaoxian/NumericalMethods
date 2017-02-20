#ifndef BinomialBarrierOption_HPP
#define BinomialBarrierOption_HPP

#include "BinomialTree.hpp"
#include "Binomial_BlackScholes.hpp"

class BinomialDownOutOption : public BinomialTree
{
protected:
	double barrier;
public:
	BinomialDownOutOption(std::tuple<double,double,double,double,double,double> data, int n_, bool isCall_, bool isEuro_,double barrier_):
		BinomialTree(data,n_,isCall_,isEuro_),barrier(barrier_){}
	double BarrierExact() {
		double a = (r - q) / sigma / sigma - 0.5;
		return BSPrice(s0, sigma, q, r, k, t, isCall) - std::pow(barrier / s0, 2.0*a)*BSPrice(barrier*barrier / s0, sigma, q, r, k, t, isCall);
	}
	virtual std::tuple<double, double, double, double> Pricing() const {
		std::tuple<double, double, double, double> result;
		std::vector<Vector> vprice, stockprice;
		vprice.resize(3);
		stockprice.resize(3);
		Vector V(n + 1);
		double dt = t / double(n);
		double u = std::exp(sigma*std::sqrt(dt));
		double d = 1.0 / u;
		double pu = (std::exp((r - q)*dt) - d) / (u - d);
		double pd = 1 - pu;

		for (int i = 0; i <= n; ++i) {
			double scurrent = s0*std::pow(u, n - i)*std::pow(d, i);
			if (isCall == true)
				V(i) = std::max(scurrent - k, 0.0);
			else
				V(i) = std::max(k - scurrent, 0.0);
			V(i) = scurrent < barrier ? 0.0 : V(i);
		}
		for (int j = n - 1; j >= 0; --j) {
			for (int i = 0; i <= j; ++i) {
				double scurrent = s0*std::pow(u, j - i)*std::pow(d, i);
				double payoff = isCall ? std::max(scurrent - k, 0.0) : std::max(k - scurrent, 0.0);
				V(i) = isEuro ? std::exp(-r*dt)*(pu*V(i) + pd*V(i + 1)) : std::max(std::exp(-r*dt)*(pu*V(i) + pd*V(i + 1)), payoff);
				V(i) = scurrent < barrier ? 0.0 : V(i);
			}
			if (j <= 2) {
				Vector v(j + 1);
				Vector stock(j + 1);
				for (int mn = 0; mn <= j; ++mn) {
					v(mn) = V(mn);
					stock(mn) = s0*std::pow(u, j - mn)*std::pow(d, mn);
				}
				vprice[j] = v;
				stockprice[j] = stock;
			}
		}
		result = std::make_tuple(V(0), Delta(vprice, stockprice), Gamma(vprice, stockprice), Theta(vprice, stockprice, dt));
		return result;
	}
	virtual std::tuple<double, double> PriceError() {
		double price = std::get<0>(Pricing());
		double error = std::abs(price - BarrierExact());
		std::tuple<double, double> result(price, error);
		return result;
	}
};

#endif
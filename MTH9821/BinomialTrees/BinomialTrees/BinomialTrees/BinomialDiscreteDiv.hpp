#ifndef BinomialDiscreteDiv_HPP
#define BinomialDiscreteDiv_HPP

#include "BinomialTree.hpp"
#include<vector>
//assume no dividend at maturity
class BinomialProportionDividend : public BinomialTree
{
protected:
	std::vector<double> DivTimes;
	std::vector<double> Dividends;
public:
	BinomialProportionDividend(std::vector<double>& DivTimes_, std::vector<double> Dividends_, std::tuple<double, double, double, double, double, double> data, int n_, bool isCall_,bool isEuro_) :
		BinomialTree(data, n_, isCall_, isEuro_),DivTimes(DivTimes_),Dividends(Dividends_) {}
	virtual std::tuple<double, double, double, double> Pricing(){
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
			double seffect = scurrent;
			for (int j = 0; j < DivTimes.size(); ++j)
				seffect = seffect*(1 - Dividends[j]);
			if (isCall == true)
				V(i) = std::max(seffect - k, 0.0);
			else
				V(i) = std::max(k - seffect, 0.0);
		}
		for (int j = n - 1; j >= 0; --j) {
			double currentTime = dt*j;
			for (int i = 0; i <= j; ++i) {
				double scurrent = s0*std::pow(u, j - i)*std::pow(d, i);
				double scurrent2 = scurrent;
				for (int k = 0; k < DivTimes.size(); ++k) {
					if (std::abs(DivTimes[k] - currentTime) > 0.00000001 && DivTimes[k] < currentTime) {
						scurrent *= (1 - Dividends[k]);
						scurrent2 *= (1-Dividends[k]);
					}
					else if (std::abs(DivTimes[k]-currentTime)<=0.00000001) {
						scurrent2 *= (1 - Dividends[k]);
					}
				}
				double payoff1 = isCall ? std::max(scurrent - k, 0.0) : std::max(k - scurrent, 0.0);
				double payoff2 = isCall ? std::max(scurrent2 - k, 0.0) : std::max(k - scurrent2, 0.0);
				double payoff = std::max(payoff1, payoff2);
				V(i) = isEuro ? std::exp(-r*dt)*(pu*V(i) + pd*V(i + 1)) : std::max(std::exp(-r*dt)*(pu*V(i) + pd*V(i + 1)), payoff);
			}
			if (j <= 2) {
				Vector v(j + 1);
				Vector stock(j + 1);
				for (int mn = 0; mn <= j; ++mn) {
					v(mn) = V(mn);
					stock(mn) = s0*std::pow(u, j - mn)*std::pow(d, mn);
					for (int k = 0; k < DivTimes.size(); ++k) {
						if (DivTimes[k] < currentTime) {
							stock(mn) *= (1 - Dividends[k]);
						}
					}
				}
				vprice[j] = v;
				stockprice[j] = stock;
			}
		}
		result = std::make_tuple(V(0), Delta(vprice, stockprice), Gamma(vprice, stockprice), Theta(vprice, stockprice, dt));
		return result;
	}
};

class BinomialFixedDividend : public BinomialTree
{
protected:
	std::vector<double> DivTimes;
	std::vector<double> Dividends;
public:
	BinomialFixedDividend(std::vector<double>& DivTimes_, std::vector<double> Dividends_, std::tuple<double, double, double, double, double, double> data, int n_, bool isCall_, bool isEuro_) :
		BinomialTree(data, n_, isCall_, isEuro_), DivTimes(DivTimes_), Dividends(Dividends_) {}
	virtual std::tuple<double, double, double, double> Pricing() {
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
		double s0NoDiv = s0;
		for (int i = 0; i < DivTimes.size(); ++i) {
			s0NoDiv -= Dividends[i] * std::exp(-r*DivTimes[i]);
		}
		for (int i = 0; i <= n; ++i) {
			double scurrent = s0NoDiv*std::pow(u, n - i)*std::pow(d, i);
			if (isCall == true)
				V(i) = std::max(scurrent - k, 0.0);
			else
				V(i) = std::max(k - scurrent, 0.0);
		}
		for (int j = n - 1; j >= 0; --j) {
			double currentTime = dt*j;
			for (int i = 0; i <= j; ++i) {
				double scurrent = s0NoDiv*std::pow(u, j - i)*std::pow(d, i);
				double scurrent2 = scurrent;
				for (int k = 0; k < DivTimes.size(); ++k) {
					if (DivTimes[k] > currentTime && std::abs(DivTimes[k] - currentTime) > 0.00000001) {
						scurrent += Dividends[k]*std::exp(-r*(DivTimes[k]-currentTime));
						scurrent2 += Dividends[k]*std::exp(-r*(DivTimes[k]-currentTime));
					}
					else if (std::abs(DivTimes[k] - currentTime) <= 0.00000001) {
						scurrent2 += Dividends[k];
					}
				}
				double payoff1 = isCall ? std::max(scurrent - k, 0.0) : std::max(k - scurrent, 0.0);
				double payoff2 = isCall ? std::max(scurrent2 - k, 0.0) : std::max(k - scurrent2, 0.0);
				double payoff = std::max(payoff1, payoff2);
				V(i) = isEuro ? std::exp(-r*dt)*(pu*V(i) + pd*V(i + 1)) : std::max(std::exp(-r*dt)*(pu*V(i) + pd*V(i + 1)), payoff);
			}
			if (j <= 2) {
				Vector v(j + 1);
				Vector stock(j + 1);
				for (int mn = 0; mn <= j; ++mn) {
					v(mn) = V(mn);
					stock(mn) = s0NoDiv*std::pow(u, j - mn)*std::pow(d, mn);
					for (int k = 0; k < DivTimes.size(); ++k) {
						if (DivTimes[k] > currentTime) {
							stock(mn) += Dividends[k] * std::exp(-r*(DivTimes[k] - currentTime));
						}
					}
				}
				vprice[j] = v;
				stockprice[j] = stock;
			}
		}
		result = std::make_tuple(V(0), Delta(vprice, stockprice), Gamma(vprice, stockprice), Theta(vprice, stockprice, dt));
		return result;
	}
};


#endif
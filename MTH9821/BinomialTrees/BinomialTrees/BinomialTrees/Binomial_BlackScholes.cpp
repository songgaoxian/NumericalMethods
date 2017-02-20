#include"Binomial_BlackScholes.hpp"
#include<random>

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


BinomialBlackScholes::BinomialBlackScholes(std::tuple<double, double, double, double, double, double> data,int n_, bool isCall_, bool isEuro_) : n(n_),
	s(std::get<0>(data)), sigma(std::get<1>(data)), q(std::get<2>(data)), r(std::get<3>(data)), k(std::get<4>(data)), t(std::get<5>(data)), isCall(isCall_), isEuro(isEuro_)
{}
std::tuple<double, double, double, double> BinomialBlackScholes::Pricing() {
	std::tuple<double, double, double, double> result;
	std::vector<Vector> vprice, stockprice;
	vprice.resize(3);
	stockprice.resize(3);
	Vector V(n);
	double dt = t / double(n);
	double u = std::exp(sigma*std::sqrt(dt));
	double d = 1.0 / u;
	double pu = (std::exp((r - q)*dt) - d) / (u - d);
	double pd = 1 - pu;

	for (int i = 0; i <= n - 1; ++i) {
		double currentPrice = s*std::pow(u, n - 1 - i)*std::pow(d, i);
		double bsPrice = BSPrice(currentPrice, sigma, q, r, k, dt, isCall);
		double payoff = isCall ? currentPrice - k : k - currentPrice; //do not need to compare with zero
		V(i) = isEuro ? bsPrice : std::max(bsPrice, payoff);
	}

	for (int j = n - 2; j >= 0; --j) {
		for (int i = 0; i <= j; ++i) {
			double currentPrice = s*std::pow(u, j - i)*std::pow(d, i);
			double payoff = isCall ? currentPrice - k : k - currentPrice;
			V(i) = isEuro ? std::exp(-r*dt)*(pu*V(i) + pd*V(i + 1)) : std::max(std::exp(-r*dt)*(pu*V(i) + pd*V(i + 1)),payoff);
		}
		if (j <= 2) {
			Vector v(j + 1);
			Vector stock(j + 1);
			for (int mn = 0; mn <= j; ++mn) {
				v(mn) = V(mn);
				stock(mn) = s*std::pow(u, j - mn)*std::pow(d, mn);
			}
			vprice[j] = v;
			stockprice[j] = stock;
		}
	}
	result = std::make_tuple(V(0), Delta(vprice, stockprice), Gamma(vprice, stockprice), Theta(vprice, stockprice, dt));
	return result;
}
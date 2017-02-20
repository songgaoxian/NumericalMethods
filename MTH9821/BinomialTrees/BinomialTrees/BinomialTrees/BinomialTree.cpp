//implement BinomialTree
#include"BinomialTree.hpp"
#include <iostream>

double Delta(std::vector<Vector> vdata, std::vector<Vector> stockdata) {
	Vector v1(vdata[1]), s1(stockdata[1]);
	return (v1(0) - v1(1)) / (s1(0) - s1(1));
}
double Gamma(std::vector<Vector> vdata, std::vector<Vector> stockdata) {
	Vector v2(vdata[2]), s2(stockdata[2]);
	return ((v2(0) - v2(1)) / (s2(0) - s2(1)) - (v2(1) - v2(2)) / (s2(1) - s2(2))) / 0.5 / (s2(0) - s2(2));
}
double Theta(std::vector<Vector> vdata, std::vector<Vector> stockdata, double dt) {
	Vector v2(vdata[2]), v0(vdata[0]);
	return (v2(1) - v0(0)) / 2.0 / dt;
}
double cdfNorm2(double x) {
	return 0.5 + 0.5*std::erf(x / std::sqrt(2.0));
}

std::tuple<double, double, double, double> BSPricingBT(double s, double sigma, double q, double r, double k, double t, bool isCall) {
	double d1 = (std::log(s / k) + (r - q + 0.5*sigma*sigma)*t) / sigma / std::sqrt(t);
	double d2 = d1 - sigma*std::sqrt(t);
	double price = isCall ? s*std::exp(-q*t)*cdfNorm2(d1) - k*std::exp(-r*t)*cdfNorm2(d2) : k*std::exp(-r*t)*cdfNorm2(-d2) - s*std::exp(-q*t)*cdfNorm2(-d1);
	double delta = isCall ? std::exp(-q*t)*cdfNorm2(d1) : -std::exp(-q*t)*cdfNorm2(-d1);
	double gamma = std::exp(-q*t) / s / sigma / std::sqrt(t) / std::sqrt(2.0 * M_PI)*std::exp(-d1*d1 / 2.0);
	double common = -s*sigma*std::exp(-q*t)*std::exp(-d1*d1 / 2.0) / 2.0 / std::sqrt(2.0*M_PI*t);
	double theta = isCall ? common + q*s*std::exp(-q*t)*cdfNorm2(d1) - r*k*std::exp(-r*t)*cdfNorm2(d2) : common - q*s*std::exp(-q*t)*cdfNorm2(-d1) + r*k*std::exp(-r*t)*cdfNorm2(-d2);
	std::tuple<double, double, double, double> result(price, delta, gamma, theta);
	return result;
}

double BSPriceBT(double s, double sigma, double q, double r, double k, double t, bool isCall) {
	double d1 = (std::log(s / k) + (r - q + 0.5*sigma*sigma)*t) / sigma / std::sqrt(t);
	double d2 = d1 - sigma*std::sqrt(t);
	return isCall ? s*std::exp(-q*t)*cdfNorm2(d1) - k*std::exp(-r*t)*cdfNorm2(d2) : k*std::exp(-r*t)*cdfNorm2(-d2) - s*std::exp(-q*t)*cdfNorm2(-d1);
}

Vector VarRedResult(Vector us_approx, Vector euro_approx, Vector euro_exact) {
	Vector VarRed = us_approx - (euro_approx - euro_exact);
	return VarRed;
}

BinomialTree::BinomialTree(std::tuple<double, double, double, double, double, double> data, int n_, bool isCall_, bool isEuro_):
	s0(std::get<0>(data)),sigma(std::get<1>(data)),q(std::get<2>(data)),r(std::get<3>(data)),k(std::get<4>(data)), t(std::get<5>(data)),n(n_),isCall(isCall_),isEuro(isEuro_){}

std::tuple<double,double,double,double> BinomialTree::Pricing() const {
	std::tuple<double,double,double,double> result;
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
	}
	//std::cout << "\n";
	for (int j = n - 1; j >= 0; --j) {
		for (int i = 0; i <= j; ++i) {
			double scurrent = s0*std::pow(u, j - i)*std::pow(d, i);
			double payoff = isCall ? std::max(scurrent - k, 0.0) : std::max(k - scurrent, 0.0);
			V(i) = isEuro ? std::exp(-r*dt)*(pu*V(i) + pd*V(i + 1)) : std::max(std::exp(-r*dt)*(pu*V(i) + pd*V(i + 1)), payoff);
		}
		//std::cout << "\n";
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
	result = std::make_tuple(V(0), Delta(vprice, stockprice), Gamma(vprice, stockprice), Theta(vprice, stockprice,dt));
	return result;
}

std::vector<Vector> BinomialTree::Pricing2(bool isBS) const {
	std::tuple<double, double, double, double> result;
	std::vector<Vector> vprice, stockprice, vexactPrice;
	vprice.resize(3);
	stockprice.resize(3);
	vexactPrice.resize(3);
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
	}
	//std::cout << "\n";
	for (int j = n - 1; j >= 0; --j) {
		for (int i = 0; i <= j; ++i) {
			double scurrent = s0*std::pow(u, j - i)*std::pow(d, i);
			double payoff = isCall ? std::max(scurrent - k, 0.0) : std::max(k - scurrent, 0.0);
			V(i) = isEuro ? std::exp(-r*dt)*(pu*V(i) + pd*V(i + 1)) : std::max(std::exp(-r*dt)*(pu*V(i) + pd*V(i + 1)), payoff);
		}
		//std::cout << "\n";
		if (j <= 2) {
			Vector v(j + 1);
			Vector stock(j + 1);
			Vector exactBS(j + 1);
			for (int mn = 0; mn <= j; ++mn) {
				v(mn) = V(mn);
				stock(mn) = s0*std::pow(u, j - mn)*std::pow(d, mn);
				double currentT = t - dt*j;
				exactBS(mn) = BSPriceBT(stock(mn), sigma, q, r, k, currentT, isCall);
			}
			vprice[j] = v;
			stockprice[j] = stock;
			vexactPrice[j] = exactBS;
		}
	}
	if (isBS) return vexactPrice;
	return vprice;
}
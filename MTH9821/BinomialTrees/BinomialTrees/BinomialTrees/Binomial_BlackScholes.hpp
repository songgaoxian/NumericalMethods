#ifndef Binomial_BlackScholes_HPP
#define Binomial_BlackScholes_HPP
#include"BinomialTree.hpp"
#include<cmath>

//calculate normal cdf
double cdfNorm(double x);
std::tuple<double, double, double, double> BSPricing(double s, double sigma, double q, double r, double k, double t, bool isCall);
double BSPrice(double s, double sigma, double q, double r, double k, double t, bool isCall);

class BinomialBlackScholes {
private:
	double s, sigma, q, r, k, t;
	int n;
	bool isCall, isEuro;
public:
	BinomialBlackScholes(std::tuple<double, double, double, double, double, double> data,int n_, bool isCall_, bool isEuro_);
	std::tuple<double, double, double, double> Pricing();
	void N(int s) { n = s; }
	void Sigma(double s) { sigma = s; }
	virtual std::tuple<double, double, double, double> VarRedPricing() {
		std::tuple<double, double, double, double, double, double> inputdata(s, sigma, q, r, k, t);
		auto exact = BSPricing(s, sigma, q, r, k, t, isCall);
		BinomialBlackScholes bt_euro(inputdata, n, isCall, true);
		auto euro_p = bt_euro.Pricing();
		auto us_p = Pricing();
		double error = std::get<0>(euro_p) - std::get<0>(exact);
		double us_p1 = std::get<0>(us_p) - error;
		error = std::get<1>(euro_p) - std::get<1>(exact);
		double us_delta = std::get<1>(us_p) - error;
		error = std::get<2>(euro_p) - std::get<2>(exact);
		double us_gamma = std::get<2>(us_p) - error;
		error = std::get<3>(euro_p) - std::get<3>(exact);
		double us_theta = std::get<3>(us_p) - error;
		std::tuple<double, double, double, double> result(us_p1, us_delta, us_gamma, us_theta);
		return result;
	}
};


#endif
#ifndef Trinomial_BlackScholes_HPP
#define Trinomial_BlackScholes_HPP
#include"TrinomialTree.hpp"
#include<cmath>

//calculate normal cdf
double cdfNorm(double x);
std::tuple<double, double, double, double> BSPricing(double s, double sigma, double q, double r, double k, double t, bool isCall);
double BSPrice(double s, double sigma, double q, double r, double k, double t, bool isCall);

class TrinomialBlackScholes {
private:
	double s, sigma, q, r, k, t;
	int n;
	bool isCall, isEuro;
public:
	TrinomialBlackScholes(std::tuple<double, double, double, double, double, double> data, int n_, bool isCall_, bool isEuro_);
	std::tuple<double, double, double, double> Pricing();
	void N(int s) { n = s; }
	void Sigma(double s) { sigma = s; }
};


#endif

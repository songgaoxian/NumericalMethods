#ifndef TrinomialBSR_HPP
#define TrinomialBSR_HPP
#include"Trinomial_BlackScholes.hpp"

class TrinomialBSR {
private:
	double s, sigma, q, r, k, t;
	int n;
	bool isCall, isEuro;
	TrinomialBlackScholes bbs;
public:
	TrinomialBSR(std::tuple<double, double, double, double, double, double> data, int n_, bool isCall_, bool isEuro_);
	std::tuple<double, double, double, double> Pricing();
	void N(int s) { n = s; bbs.N(s); }
	void Sigma(double s) { sigma = s; bbs.Sigma(s); }
};

#endif

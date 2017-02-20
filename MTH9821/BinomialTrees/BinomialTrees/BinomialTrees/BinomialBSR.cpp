#include"BinomialBSR.hpp"

BinomialBSR::BinomialBSR(std::tuple<double, double, double, double, double, double> data, int n_, bool isCall_, bool isEuro_) : n(n_),bbs(data,n_,isCall_,isEuro_),
s(std::get<0>(data)), sigma(std::get<1>(data)), q(std::get<2>(data)), r(std::get<3>(data)), k(std::get<4>(data)), t(std::get<5>(data)), isCall(isCall_), isEuro(isEuro_)
{}

std::tuple<double, double, double, double> BinomialBSR::Pricing() {
	auto result1 = bbs.Pricing();
	bbs.N(n / 2);
	auto result2 = bbs.Pricing();
	double price = 2 * std::get<0>(result1) - std::get<0>(result2);
	double delta = 2 * std::get<1>(result1) - std::get<1>(result2);
	double gamma = 2 * std::get<2>(result1) - std::get<2>(result2);
	double theta = 2 * std::get<3>(result1) - std::get<3>(result2);
	std::tuple<double, double, double, double> result(price, delta, gamma, theta);
	return result;
}
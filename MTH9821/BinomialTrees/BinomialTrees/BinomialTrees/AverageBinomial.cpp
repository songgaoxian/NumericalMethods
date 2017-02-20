#include"AverageBinomial.hpp"

AverageBinomial::AverageBinomial(std::tuple<double, double, double, double, double, double> data, int n_, bool isCall_, bool isEuro_) :
	bt(data, n_,isCall_, isEuro_) {}

std::tuple<double,double,double,double> AverageBinomial::Pricing() {
	auto result1 = bt.Pricing();
	int currentN = bt.N();
	bt.N(currentN + 1);
	auto result2 = bt.Pricing();
	std::tuple<double, double, double, double> result(0.5*(std::get<0>(result1) + std::get<0>(result2)), 0.5*(std::get<1>(result1) + std::get<1>(result2)), 0.5*(std::get<2>(result1) + std::get<2>(result2)), 0.5*(std::get<3>(result1) + std::get<3>(result2)));
	return result;
}
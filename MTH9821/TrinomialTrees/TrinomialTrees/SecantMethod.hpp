#ifndef SecantMethod_HPP
#define SecantMethod_HPP

#include"TrinomialTree.hpp"
#include"Trinomial_BlackScholes.hpp"
#include"TrinomialBSR.hpp"
#include"TrinomialBarrierOption.hpp"
#include<cmath>
#include<iostream>

class SecantVolatility
{
private:
	TrinomialTree bt;
	double sigma0, sigma1;
	double mkt_price;
	double tol;
public:
	SecantVolatility(TrinomialTree bt_, double sigma0_, double sigma1_, double mkt_price_, double tol_) :
		bt(bt_), sigma0(sigma0_), sigma1(sigma1_), mkt_price(mkt_price_), tol(tol_) {}
	double ImpliedVolAbsolute() {
		double x0 = sigma0;
		double x1 = sigma1;
		bt.Sigma(x0);
		double f0 = std::get<0>(bt.Pricing()) - mkt_price;
		bt.Sigma(x1);
		double f1 = std::get<0>(bt.Pricing()) - mkt_price;
		int n = 0;
		while (std::abs(f1) >= tol) {
			double x = x1 - f1*(x1 - x0) / (f1 - f0);
			bt.Sigma(x);
			double f = std::get<0>(bt.Pricing()) - mkt_price;
			f0 = f1;
			f1 = f;
			x0 = x1;
			x1 = x;
			++n;
		}
		std::cout << "runs:" << n << "\n";
		return x1;
	}
	double ImpliedVolConsecutive() {
		double x0 = sigma0;
		double x1 = sigma1;
		bt.Sigma(x0);
		double f0 = std::get<0>(bt.Pricing()) - mkt_price;
		bt.Sigma(x1);
		double f1 = std::get<0>(bt.Pricing()) - mkt_price;
		int n = 0;
		while (std::abs(x1 - x0) >= tol) {
			double x = x1 - f1*(x1 - x0) / (f1 - f0);
			bt.Sigma(x);
			double f = std::get<0>(bt.Pricing()) - mkt_price;
			f0 = f1;
			f1 = f;
			x0 = x1;
			x1 = x;
			++n;
		}
		std::cout << "runs:" << n << "\n";
		return x1;
	}
};


#endif
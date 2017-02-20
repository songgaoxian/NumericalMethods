#ifndef TrinomialBarrierOption_HPP
#define TrinomialBarrierOption_HPP

#include "TrinomialBSR.hpp"

class TrinomialDownOut : public TrinomialTree
{
protected:
	double Bup, Blow, pay;
public:
	TrinomialDownOut(double Bup_,double Blow_,double pay_,std::tuple<double, double, double, double, double, double> data, int n_, bool isCall_, bool isEuro_):
		TrinomialTree(data,n_,isCall_,isEuro_),Bup(Bup_),Blow(Blow_),pay(pay_){}
	virtual std::tuple<double, double, double, double> Pricing() const{
		std::tuple<double, double, double, double> result;
		std::vector<Vector> vprice, stockprice;
		vprice.resize(3);
		stockprice.resize(3);
		Vector V(n * 2 + 1);
		double dt = t / double(n);
		double u = std::exp(sigma*std::sqrt(3.0*dt));
		double d = 1.0 / u;
		double pu = 1.0 / 6.0 + (r - q - sigma*sigma*0.5)*std::sqrt(dt / 12.0 / sigma / sigma);
		double pm = 2.0 / 3.0;
		double pd = 1 - pu - pm;

		for (int i = 0; i <= n * 2; ++i) {
			double scurrent = s0*std::pow(u, n - i);
			V(i) = isCall ? std::max(scurrent - k, 0.0) : std::max(k - scurrent, 0.0);
			if (scurrent >= Bup || scurrent <= Blow) { V(i) = pay; }
		}
		for (int j = n - 1; j >= 0; --j) {
			for (int i = 0; i <= j * 2; ++i) {
				double scurrent = s0*std::pow(u, j - i);
				double EuroValue = std::exp(-r*dt)*(pu*V(i) + pm*V(i + 1) + pd*V(i + 2));
				double ExerciseVal = isCall ? std::max(scurrent - k, 0.0) : std::max(k - scurrent, 0.0);
				V(i) = isEuro ? EuroValue : std::max(EuroValue, ExerciseVal);
				if (scurrent >= Bup || scurrent <= Blow) { V(i) = pay; }
			}
			if (j <= 2) {
				Vector v(j * 2 + 1);
				Vector stock(j * 2 + 1);
				for (int mn = 0; mn <= j * 2; ++mn) {
					v(mn) = V(mn);
					stock(mn) = s0*std::pow(u, j - mn);
				}
				vprice[j] = v;
				stockprice[j] = stock;
			}
			if (j % 1000 == 0) std::cout << j << "\n";
		}
		result = std::make_tuple(V(0), Delta(vprice, stockprice), Gamma(vprice, stockprice), Theta(vprice, stockprice, dt));
		return result;
	}
	double ExactP() {
		double barrier = Bup;
		double a = (r - q) / sigma / sigma - 0.5;
		return BSPrice(s0, sigma, q, r, k, t, isCall) - std::pow(barrier / s0, 2.0*a)*BSPrice(barrier*barrier / s0, sigma,q,r,k,t,isCall);
	}
	std::tuple<double, double> PriceError() {
		auto temp = Pricing();
		double price = std::get<0>(temp);
		double error = std::abs(price - ExactP());
		std::tuple<double, double> result(price, error);
		return result;
	}
};

#endif
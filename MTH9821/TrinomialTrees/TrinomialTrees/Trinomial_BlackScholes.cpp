#include"Trinomial_BlackScholes.hpp"
#include<random>


TrinomialBlackScholes::TrinomialBlackScholes(std::tuple<double, double, double, double, double, double> data, int n_, bool isCall_, bool isEuro_) : n(n_),
s(std::get<0>(data)), sigma(std::get<1>(data)), q(std::get<2>(data)), r(std::get<3>(data)), k(std::get<4>(data)), t(std::get<5>(data)), isCall(isCall_), isEuro(isEuro_)
{}
std::tuple<double, double, double, double> TrinomialBlackScholes::Pricing() {
	std::tuple<double, double, double, double> result;
	std::vector<Vector> vprice, stockprice;
	vprice.resize(3);
	stockprice.resize(3);
	Vector V(n*2-1);
	double dt = t / double(n);
	double u = std::exp(sigma*std::sqrt(3.0*dt));
	double d = 1.0 / u;
	double pu = 1.0 / 6.0 + (r - q - sigma*sigma*0.5)*std::sqrt(dt / 12.0 / sigma / sigma);
	double pm = 2.0 / 3.0;
	double pd = 1 - pu - pm;

	for (int i = 0; i <= n*2 - 2; ++i) {
		double currentPrice = s*std::pow(u, n - 1 - i);
		double bsPrice = BSPrice(currentPrice, sigma, q, r, k, dt, isCall);
		double payoff = isCall ? currentPrice - k : k - currentPrice; //do not need to compare with zero
		V(i) = isEuro ? bsPrice : std::max(bsPrice, payoff);
	}

	for (int j = n - 2; j >= 0; --j) {
		for (int i = 0; i <= j*2; ++i) {
			double currentPrice = s*std::pow(u, j - i);
			double payoff = isCall ? currentPrice - k : k - currentPrice;
			V(i) = isEuro ? std::exp(-r*dt)*(pu*V(i) + pm*V(i + 1)+pd*V(i+2)) : std::max(std::exp(-r*dt)*(pu*V(i) + pm*V(i + 1)+pd*V(i+2)), payoff);
		}
		if (j <= 2) {
			Vector v(j*2 + 1);
			Vector stock(j*2 + 1);
			for (int mn = 0; mn <= j*2; ++mn) {
				v(mn) = V(mn);
				stock(mn) = s*std::pow(u, j - mn);
			}
			vprice[j] = v;
			stockprice[j] = stock;
		}
	}
	result = std::make_tuple(V(0), Delta(vprice, stockprice), Gamma(vprice, stockprice), Theta(vprice, stockprice, dt));
	return result;
}
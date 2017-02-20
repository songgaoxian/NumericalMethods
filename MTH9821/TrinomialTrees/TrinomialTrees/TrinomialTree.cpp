//implement TrinomialTree
#include"TrinomialTree.hpp"

double Delta(std::vector<Vector> vdata, std::vector<Vector> stockdata) {
	Vector v1(vdata[1]), s1(stockdata[1]);
	return (v1(0) - v1(2)) / (s1(0) - s1(2));
}
double Gamma(std::vector<Vector> vdata, std::vector<Vector> stockdata) {
	Vector v2(vdata[2]), s2(stockdata[2]), s1(stockdata[1]);
	double up = (v2(0) - v2(2)) / (s2(0) - s2(2)) - (v2(2) - v2(4)) / (s2(2) - s2(4));
	return up / (s1(0) - s1(2));
}
double Theta(std::vector<Vector> vdata, std::vector<Vector> stockdata, double dt) {
	Vector v1(vdata[1]),v0(vdata[0]);
	return (v1(1) - v0(0)) / dt;
}

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

Vector VarRedResult(Vector us_approx, Vector euro_approx, Vector euro_exact) {
	Vector VarRed = us_approx - (euro_approx - euro_exact);
	return VarRed;
}

TrinomialTree::TrinomialTree(std::tuple<double, double, double, double, double, double> data, int n_, bool isCall_, bool isEuro_) :
	s0(std::get<0>(data)), sigma(std::get<1>(data)), q(std::get<2>(data)), r(std::get<3>(data)), k(std::get<4>(data)), t(std::get<5>(data)), n(n_), isCall(isCall_), isEuro(isEuro_) {}

std::tuple<double, double, double, double> TrinomialTree::Pricing() const {
	std::tuple<double, double, double, double> result;
	std::vector<Vector> vprice, stockprice;
	vprice.resize(3);
	stockprice.resize(3);
	Vector V(n*2 + 1);
	double dt = t / double(n);
	double u = std::exp(sigma*std::sqrt(3.0*dt));
	double d = 1.0 / u;
	double pu = 1.0/6.0+(r-q-sigma*sigma*0.5)*std::sqrt(dt/12.0/sigma/sigma);
	double pm = 2.0 / 3.0;
	double pd = 1 - pu - pm;

	for (int i = 0; i <= n*2; ++i) V(i) = isCall? std::max(s0*std::pow(u,n-i)-k,0.0):std::max(k-s0*std::pow(u,n-i),0.0);

	for (int j = n - 1; j >= 0; --j) {
		for (int i = 0; i <= j*2; ++i) {
			double EuroValue = std::exp(-r*dt)*(pu*V(i) + pm*V(i + 1) + pd*V(i + 2));
			double ExerciseVal = isCall ? std::max(s0*std::pow(u, j - i) - k, 0.0) : std::max(k - s0*std::pow(u, j - i), 0.0);
			V(i) = isEuro ? EuroValue : std::max(EuroValue, ExerciseVal);
		}
		if (j <= 2) {
			Vector v(j*2 + 1);
			Vector stock(j*2 + 1);
			for (int mn = 0; mn <= j*2; ++mn) {
				v(mn) = V(mn);
				stock(mn) = s0*std::pow(u, j - mn);
			}
			vprice[j] = v;
			stockprice[j] = stock;
		}
	}
	result = std::make_tuple(V(0), Delta(vprice, stockprice), Gamma(vprice, stockprice), Theta(vprice, stockprice, dt));
	return result;
}
std::vector<Vector> TrinomialTree::Pricing2(bool isBS) const {
	std::tuple<double, double, double, double> result;
	std::vector<Vector> vprice, stockprice, vexact;
	vprice.resize(3);
	stockprice.resize(3);
	vexact.resize(3);
	Vector V(n * 2 + 1);
	double dt = t / double(n);
	double u = std::exp(sigma*std::sqrt(3.0*dt));
	double d = 1.0 / u;
	double pu = 1.0 / 6.0 + (r - q - sigma*sigma*0.5)*std::sqrt(dt / 12.0 / sigma / sigma);
	double pm = 2.0 / 3.0;
	double pd = 1 - pu - pm;

	for (int i = 0; i <= n * 2; ++i) V(i) = isCall ? std::max(s0*std::pow(u, n - i) - k, 0.0) : std::max(k - s0*std::pow(u, n - i), 0.0);

	for (int j = n - 1; j >= 0; --j) {
		for (int i = 0; i <= j * 2; ++i) {
			double EuroValue = std::exp(-r*dt)*(pu*V(i) + pm*V(i + 1) + pd*V(i + 2));
			double ExerciseVal = isCall ? std::max(s0*std::pow(u, j - i) - k, 0.0) : std::max(k - s0*std::pow(u, j - i), 0.0);
			V(i) = isEuro ? EuroValue : std::max(EuroValue, ExerciseVal);
		}
		if (j <= 2) {
			Vector v(j * 2 + 1);
			Vector stock(j * 2 + 1);
			Vector vBs(j * 2 + 1);
			for (int mn = 0; mn <= j * 2; ++mn) {
				v(mn) = V(mn);
				stock(mn) = s0*std::pow(u, j - mn);
				double currentT = t - dt*j;
				vBs(mn) = BSPrice(stock(mn), sigma, q, r, k, currentT, isCall);
			}
			vprice[j] = v;
			stockprice[j] = stock;
			vexact[j] = vBs;
		}
	}
	if (isBS) return vexact;
	return vprice;
}
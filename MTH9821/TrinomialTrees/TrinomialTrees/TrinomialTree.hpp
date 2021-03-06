//implement binomial tree pricing
#ifndef TrinomialTree_HPP
#define TrinomialTree_HPP
#include<Eigen/Dense>
#include<vector>
#include<tuple>
#include<cmath>

using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;

double Delta(std::vector<Vector> vdata, std::vector<Vector> stockdata);
double Gamma(std::vector<Vector> vdata, std::vector<Vector> stockdata);
double Theta(std::vector<Vector> vdata, std::vector<Vector> stockdata, double dt);
Vector VarRedResult(Vector us_approx, Vector euro_approx, Vector euro_exact);

//calculate normal cdf
double cdfNorm(double x);
std::tuple<double, double, double, double> BSPricing(double s, double sigma, double q, double r, double k, double t, bool isCall);
double BSPrice(double s, double sigma, double q, double r, double k, double t, bool isCall);


class TrinomialTree {
protected:
	double s0, sigma, q, r, k, t;
	int n;
	bool isEuro, isCall;
public:
	TrinomialTree(std::tuple<double, double, double, double, double, double> data, int n_, bool isCall_, bool isEuro_);
	void N(int s) { n = s; }
	double DT() const { return t / double(n); }
	//calculate binomial pricing
	virtual std::tuple<double, double, double, double> Pricing() const;
	void Sigma(double s) { sigma = s; }
	virtual std::vector<Vector> Pricing2(bool isBS) const;

	int StepsForConverge(double tol) {
		double p1 = std::get<0>(Pricing());
		n *= 2;
		double p2 = std::get<0>(Pricing());
		while (std::abs(p1 - p2) >= tol) {
			p1 = p2;
			n *= 2;
			p2 = std::get<0>(Pricing());
		}
		return n;
	}

	int StepsForConvergeUSVarRed(double tol) {
		std::tuple<double, double, double, double, double, double> inputdata(s0, sigma, q, r, k, t);
		double euprice = BSPrice(s0, sigma, q, r, k, t, isCall);
		TrinomialTree bt_euro(inputdata, n, isCall, true);
		double p1 = std::get<0>(Pricing()) - (std::get<0>(bt_euro.Pricing()) - euprice);
		n *= 2;
		bt_euro.N(n);
		double p2 = std::get<0>(Pricing()) - (std::get<0>(bt_euro.Pricing()) - euprice);
		while (std::abs(p1 - p2) >= tol) {
			p1 = p2;
			n *= 2;
			bt_euro.N(n);
			p2 = std::get<0>(Pricing()) - (std::get<0>(bt_euro.Pricing()) - euprice);
		}
		return n;
	}
	virtual std::tuple<double, double, double, double> VarRedPricing() {
		std::tuple<double, double, double, double, double, double> inputdata(s0, sigma, q, r, k, t);
		auto exact = BSPricing(s0, sigma, q, r, k, t, isCall);
		TrinomialTree bt_euro(inputdata, n, isCall, true);
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

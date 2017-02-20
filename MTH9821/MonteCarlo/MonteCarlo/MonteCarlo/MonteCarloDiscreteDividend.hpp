#ifndef MonteCarloDiscreteDividend_HPP
#define MonteCarloDiscreteDividend_HPP

#include "MonteCarlo.hpp"
#include<functional>
#include<vector>

using DividendFunct = std::function<double(double)>;
double Average(std::vector<double>& input) {
	double sum = 0;
	for (auto it = input.begin(); it != input.end(); ++it)
		sum += *it;
	return sum / double(input.size());
}
std::vector<double> ControlVariate(std::vector<double>& prime, std::vector<double>& control_v, double BS) {
	double primeAverage = Average(prime);
	double controlAverage = Average(control_v);
	double b_up = 0, b_down = 0;
	for (int i = 0; i < prime.size(); ++i) {
		b_up += (prime[i] - primeAverage)*(control_v[i] - controlAverage);
		b_down += (control_v[i] - controlAverage)*(control_v[i] - controlAverage);
	}
	double b = b_up / b_down;
	std::vector<double> result;
	for (int i = 0; i < prime.size(); ++i) {
		result.push_back(prime[i] - b*(control_v[i] - BS));
	}
	return result;
}

class MonteCarloDiscreteDividend : public MonteCarlo
{
protected:
	std::vector<double> time_dividends;
	std::vector<DividendFunct> dividendFunctions;
	bool isCall;
public:
	MonteCarloDiscreteDividend(long size_, std::tuple<double, double, double, double, double, double>& stck, int type, std::vector<double> times_, std::vector<DividendFunct> functs_,bool isCall_) :
		MonteCarlo(size_*(time_dividends.size()+1)*6, stck, type), time_dividends(times_), dividendFunctions(functs_), isCall(isCall_) {
		std::cout << "num of normals=" << normal.size() << "\n";
	}
	std::tuple<double, double> PriceDelta(long n) {
		long count = 0;
		std::vector<double> P, Delta;
		std::vector<double> Pbar, Deltabar;
		double s0 = std::get<0>(stock);
		double k = std::get<1>(stock);
		double t = std::get<2>(stock);
		double sigma = std::get<3>(stock);
		double q = std::get<4>(stock);
		double r = std::get<5>(stock);
		//calculate delta
		std::function<double(double, double)> deltafunct = [=](double st, double stbar) {
			double ind = k > st ? -1.0 : 0.0;
			return ind*std::exp(-r*t)*stbar / s0*(1 - 0.01);
		};
		std::function<double(double)> deltafunctBar = [=](double stbar) {
			double ind = k > stbar ? -1.0 : 0.0;
			return ind*std::exp(-r*t)*stbar / s0;
		};
		std::vector<double> zs(time_dividends.size()+1);
		for (int i = 0; i < n; ++i) {
			double siDiv = s0;
			double siNoDiv = s0;
			bool valid = true;
			for (int j = 0; j < time_dividends.size(); ++j) {
				double t1, t2;
				if (j == 0) { t1 = 0; }
				else { t1 = time_dividends[j - 1]; }
				t2 = time_dividends[j];
				zs[j] = normal[count++];
				double sPrevDiv = siDiv*std::exp((r - sigma*sigma*0.5)*(t2 - t1) + sigma*std::sqrt(t2 - t1)*zs[j]);
				siDiv = dividendFunctions[j](sPrevDiv);
				if (siDiv < 0) valid = false;
				siNoDiv = siNoDiv*std::exp((r - sigma*sigma*0.5)*(t2 - t1) + sigma*std::sqrt(t2 - t1)*zs[j]);
			}
			if (valid) {
				double t1 = time_dividends.back();
				zs[time_dividends.size()] = normal[count++];
				siDiv = siDiv*std::exp((r - sigma*sigma*0.5)*(t - t1) + sigma*std::sqrt(t - t1)*zs[time_dividends.size()]);
				siNoDiv = siNoDiv*std::exp((r - sigma*sigma*0.5)*(t - t1) + sigma*std::sqrt(t - t1)*zs[time_dividends.size()]);
				double payoff = isCall ? std::max(siDiv - k, 0.0) : std::max(k - siDiv, 0.0);
				double payoffBar = isCall ? std::max(siNoDiv - k, 0.0) : std::max(k - siNoDiv, 0.0);
				P.push_back(payoff);
				Pbar.push_back(payoffBar);
				Delta.push_back(deltafunct(siDiv, siNoDiv));
				Deltabar.push_back(deltafunctBar(siNoDiv));
			}
		}
		auto temp = BSPricing();
		double bs_price = isCall ? std::get<0>(temp) : std::get<3>(temp);
		double bs_delta = isCall ? std::get<1>(temp) : std::get<4>(temp);
		auto Wi = ControlVariate(P, Pbar, bs_price);
		auto Deltas = ControlVariate(Delta, Deltabar, bs_delta);
		std::tuple<double, double> result(Average(Wi), Average(Deltas));
		return result;
	}		
};   



#endif
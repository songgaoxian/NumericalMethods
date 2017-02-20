#ifndef MonteCarlo_HPP
#define MonteCarlo_HPP

#include<tuple>
#include<algorithm>
#include<iostream>
#include "Random.hpp"
#include<random>

#define PI 3.14159265358979323846

class MonteCarlo {
protected:
	long max_size;
	RandomGenerator rg;
	std::vector<double> normal;
	std::tuple<double, double, double, double, double, double> stock;
public:
	MonteCarlo(long size_,std::tuple<double,double,double,double,double,double>& stck, int type):rg(size_),max_size(size_),stock(stck){
		if (type == 0) {
			normal = rg.InverseTransform();
		}
		else if (type == 1) {
			normal = rg.AcceptReject();
		}
		else {
			normal = rg.Box_Muller();
		}
	}
	std::tuple<double, double, double, double, double, double> BSPricing() {
		double s0 = std::get<0>(stock);
		double k = std::get<1>(stock);
		double t = std::get<2>(stock);
		double sigma = std::get<3>(stock);
		double q = std::get<4>(stock);
		double r = std::get<5>(stock);
		double d1 = (std::log(s0 / k) + (r - q + 0.5*sigma*sigma)*t) / (sigma*std::sqrt(t));
		double d2 = d1 - sigma*std::sqrt(t);
		double c = s0*std::exp(-q*t)*cdfNorm(d1) - k*std::exp(-r*t)*cdfNorm(d2);
		double p = k*std::exp(-r*t)*cdfNorm(-d2) - s0*std::exp(-q*t)*cdfNorm(-d1);
		double deltac = std::exp(-q*t)*cdfNorm(d1);
		double deltap = deltac - std::exp(-q*t);
		double vegac = std::sqrt(t)*s0*std::exp(-q*t)*1.0 / std::sqrt(2.0*PI)*std::exp(-d1*d1 / 2.0);
		double vegap = vegac;
		std::tuple<double, double, double, double, double, double> result(c, deltac, vegac, p, deltap, vegap);
		return result;
	}

	std::tuple<double, double, double, double, double, double> CalculatePlain(long N) {
		//another version
		//N = normal.size();
		std::cout << "normal samples=" << N << "\n";
		double s, c=0, deltac=0, vegac=0, p=0, deltap=0, vegap=0;
		double s0 = std::get<0>(stock);
		double k = std::get<1>(stock);
		double t = std::get<2>(stock);
		double sigma = std::get<3>(stock);
		double q = std::get<4>(stock);
		double r = std::get<5>(stock);
		for (long i = 0; i < N; ++i) {
			s = s0*std::exp((r - q - 0.5*sigma*sigma)*t + sigma*std::sqrt(t)*normal[i]);
			c += std::exp(-r*t)*std::max(0.0, s - k) / double(N);
			p += std::exp(-r*t)*std::max(0.0, k - s) / double(N);
			if (s > k) {
				deltac += std::exp(-r*t)*s / s0 / double(N);
				vegac += s*std::exp(-r*t)*(-sigma*t + std::sqrt(t)*normal[i]) / double(N);
				deltap += 0.0 / double(N);
				vegap += 0.0 / double(N);
			}
			else {
				deltac += 0.0 / double(N);
				vegac += 0.0 / double(N);
				deltap += -std::exp(-r*t)*s / s0 / double(N);
				vegap += -s*std::exp(-r*t)*(-sigma*t + std::sqrt(t)*normal[i]) / double(N);
			}
		}
		std::tuple<double, double, double, double, double, double> result(c, deltac, vegac, p, deltap, vegap);
		return result;
	}

	double BSCall(double S0) {
		double s0 = S0;
		double k = std::get<1>(stock);
		double t = std::get<2>(stock);
		double sigma = std::get<3>(stock);
		double q = std::get<4>(stock);
		double r = std::get<5>(stock);
		double d1 = (std::log(s0 / k) + (r - q + 0.5*sigma*sigma)*t) / (sigma*std::sqrt(t));
		double d2 = d1 - sigma*std::sqrt(t);
		double c = s0*std::exp(-q*t)*cdfNorm(d1) - k*std::exp(-r*t)*cdfNorm(d2);
		//double p=k*std::exp(-r*t)*cdfNorm(-d2)-s0*std::exp(-q*t)*cdfNorm(-d1);
		return c;
	}
	double DownAndOutCall(double B) {
		double s0 = std::get<0>(stock);
		double k = std::get<1>(stock);
		double t = std::get<2>(stock);
		double sigma = std::get<3>(stock);
		double q = std::get<4>(stock);
		double r = std::get<5>(stock);
		double result = BSCall(s0) - std::pow(s0 / B, 1 - 2 * (r - q) / (sigma*sigma))*BSCall(B*B / s0);
		return result;
	}
	//N is total random normal samples
	std::tuple<double,double,double,double> CalculateDownAndOutCall(double B,long m, long N) {
		long n = std::floor(N/ m);
		double v = 0;
		double s;
		long j = 0;
		double s0 = std::get<0>(stock);
		double k = std::get<1>(stock);
		double t = std::get<2>(stock);
		double sigma = std::get<3>(stock);
		double q = std::get<4>(stock);
		double r = std::get<5>(stock);
		//another version
		//m = std::ceil(std::pow(N, 1.0 / 3.0)*std::pow(t, 2.0 / 3.0));
		//n = std::floor(N / m);
		double dt = t / double(m);
		long count = 0;
		for (long i = 0; i <n; ++i) {
			j = 0;
			s = s0;
			bool isValid = true;
			while (j < m) {
				++j;
				s = s*std::exp((r-q- 0.5*sigma*sigma)*dt + sigma*std::sqrt(dt)*normal[count++]);
				if (s <= B) { isValid = false; }
			}
			if (isValid) { v += std::max(s - k, 0.0)*std::exp(-r*t) / double(n); }
		}
		std::tuple<double, double, double,double> result(m, n, v, std::abs(v - DownAndOutCall(B)));
		return result;
	}

};


#endif
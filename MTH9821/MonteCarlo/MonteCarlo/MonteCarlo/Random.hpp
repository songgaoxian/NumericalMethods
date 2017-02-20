#ifndef Random_HPP
#define Random_HPP

#define _USE_MATH_DEFINES
#include<vector>
#include<climits>
#include<cmath>

double ApproximateInverse(double u) {
	double y = u - 0.5, r, x;
	double a0 = 2.50662823884, a1 = -18.61500062529, a2 = 41.39119773534, a3 = -25.44106049637;
	double b0 = -8.4735109309, b1 = 23.08336743743, b2 = -21.06224101826, b3 = 3.13082909833;
	double c0 = 0.3374754822726147, c1 = 0.9761690190917186, c2 = 0.1607979714918209, c3 = 0.0276438810333863;
	double c4 = 0.0038405729373609, c5 = 0.0003951896511919, c6 = 0.0000321767881768;
	double c7 = 0.0000002888167364, c8 = 0.0000003960315187;
	if (std::abs(y) < 0.42) {
		r = y*y;
		x = y*(((a3*r + a2)*r + a1)*r + a0) / ((((b3*r + b2)*r + b1)*r + b0)*r + 1);
	}
	else {
		r = u;
		if (y > 0) r = 1 - u;
		r = std::log(-std::log(r));
		x = c0 + r*(c1 + r*(c2 + r*(c3 + r*(c4 + r*(c5 + r*(c6 + r*(c7 + r*c8)))))));
		if (y < 0) x = -x;
	}
	return x;
}

double cdfNorm(double x) {
	return 0.5 + 0.5*std::erf(x / std::sqrt(2.0));
}

class RandomGenerator {
protected:
	std::vector<double> uniform;
	std::vector<double> normal;
	long N;
public:
	RandomGenerator(long srcN){
		N = srcN;
		long long k = std::pow(2,31)-1;
		long long a1 = 8192 * 4;
		long long a2 = 6605;
		long long x0 = 1;
		//uniform.push_back(double(x0) / double(k));
		for (long i = 1; i <= N; ++i) {
			long long temp1 = (a1*x0) % k;
			long long temp2 = (a2*x0) % k;
			x0 = (temp1 + temp2) % k;
			uniform.push_back(double(x0) / double(k));
		}
	}

	std::vector<double> GetUniform() { return uniform; }
	long GetN() { return N; }

	std::vector<double> InverseTransform() { 
		normal.clear();
		for (auto it = uniform.begin(); it != uniform.end(); ++it) normal.push_back(ApproximateInverse(*it));
		return normal;
	}

	std::vector<double> AcceptReject() {
		normal.clear();
		auto it = uniform.begin();
		double x;
		while (it != uniform.end()) {
			double u1 = *it; ++it;
			if (it == uniform.end()) return normal;
			double u2 = *it; ++it;
			if (it == uniform.end()) return normal;
			double u3 = *it; ++it;
			x = -std::log(u1);
			if (u2 <= std::exp(-0.5*(x - 1)*(x - 1))) {
				if (u3 <= 0.5) x = -x;
				normal.push_back(x);
			}
		}
		return normal;
	}
	
	std::vector<double> Box_Muller() {
		normal.clear();
		double x, u1, u2, y, z1, z2;
		auto it = uniform.begin();
		while (it != uniform.end()) {
			do {
				if (it == uniform.end()) return normal;
				u1 = *it;
				++it;
				if (it == uniform.end()) return normal;
				u2 = *it;
				u1 = 2 * u1 - 1; u2 = 2 * u2 - 1;
				x = u1*u1 + u2*u2;
				++it;
			} while (x > 1);
			y = std::sqrt(-2 * std::log(x) / x);
			z1 = u1*y; 
			z2 = u2*y;
			normal.push_back(z1);
			normal.push_back(z2);
		}
		return normal;
	}
};



#endif